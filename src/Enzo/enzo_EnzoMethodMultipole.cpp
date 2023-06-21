// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMultipole.cpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     July 19, 2022
/// @brief    Compute multipoles and pass multipoles up octree

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodMultipole::EnzoMethodMultipole (double timeStep, double theta)
  : Method(),
    timeStep_(timeStep),
    theta_(theta),
    is_volume_(-1),
    block_volume_(),
    max_volume_(0)
{ 

  cello::define_field ("density");
  cello::define_field ("acceleration_x");
  cello::define_field ("acceleration_y");
  cello::define_field ("acceleration_z");

  // Initialize default Refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("density");
  refresh->add_field("acceleration_x");
  refresh->add_field("acceleration_y");
  refresh->add_field("acceleration_z");

  // Initialize Sync scalars
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  i_sync_restrict_ = scalar_descr_sync->new_value("multipole:restrict");

  // Initialize scalars for receiving messages
  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  i_msg_prolong_ = scalar_descr_void->new_value("multipole:msg_prolong");
  for (int ic=0; ic<cello::num_children(); ic++) {
    i_msg_restrict_[ic] = scalar_descr_void->new_value("multipole:msg_restrict");
  }

  // Initialize mass, COM, and quadrupole scalars
  ScalarDescr * scalar_descr_double = cello::scalar_descr_double();
  i_mass_ = scalar_descr_double->new_value("multipole:mass");
  i_com_ = scalar_descr_double->new_value("multipole:com", 3);
  i_quadrupole_ = scalar_descr_double->new_value("multipole:quadrupole", 9);

  // Initialize Taylor coefficient scalars
  i_c1_ = scalar_descr_double->new_value("multipole:c1", 3);
  i_c2_ = scalar_descr_double->new_value("multipole:c2", 9);
  i_c3_ = scalar_descr_double->new_value("multipole:c3", 27);


  // Declare long long Block Scalar for volume and save scalar index
  is_volume_ = cello::scalar_descr_long_long()->new_value("solver_fmm_volume");


}

//----------------------------------------------------------------------

void EnzoMethodMultipole::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | timeStep_;
  p | i_sync_restrict_;
  p | i_msg_prolong_;
  PUParray(p, i_msg_restrict_, 8);

  p | i_mass_;
  p | i_com_;
  p | i_quadrupole_;
  p | i_c1_;
  p | i_c2_;
  p | i_c3_;

  p | theta_;
  p | is_volume_;
  p | block_volume_;
  p | max_volume_;
  
}

//----------------------------------------------------------------------

void EnzoMethodMultipole::compute ( Block * block) throw()
{

  Sync * sync_restrict = psync_restrict(block);
  //Sync * sync_prolong = psync_prolong(block);
  sync_restrict->set_stop(1 + cello::num_children());
  //sync_prolong->set_stop(1 + 1);
  
  int level = block->level();

  double * mass = pmass(block);
  double * com = pcom(block);
  double * quadrupole = pquadrupole(block);
  double * c1 = pc1(block);
  double * c2 = pc2(block);
  double * c3 = pc3(block);
  
  *mass = 0;

  for (int i = 0; i < 3; i++){
    com[i]  = 0;   
  }
  
  for (int i = 0; i < 9; i++){
    quadrupole[i]  = 0;
  }

  for (int i = 0; i < 3; i++){
      c1[i] = 0;
  }
  
  for (int i = 0; i < 9; i++){
      c2[i] = 0;
    }

  for (int i = 0; i < 27; i++){
      c3[i] = 0;
  }


  // if (level == 0) {    // for testing purposes

  //   for (int i = 0; i < 3; i++){
  //     c1[i] = 1;
  //   }

  //   for (int i = 0; i < 9; i++){
  //     c2[i] = 1;
  //   }

  //   c3[0] = 1; 
  // }


  // Allocate and initialize block_volume_[level - min_level] to store
  // weighted volume of blocks in different levels relative to the finest
  // (for the time being, I've set min_level = 0)
  const int num_children = cello::num_children();
  const int max_level_ = cello::hierarchy()->max_level();

  block_volume_.resize(max_level_ + 1);
  block_volume_[max_level_] = 1;

  for (int i = max_level_ - 1; i >= 0; i--) {
    block_volume_[i] = block_volume_[i+1] * num_children;
  }

  // Compute volume of the domain relative to finest-level blocks
  int nb3[3] = {1,1,1};
  cello::hierarchy()->root_blocks(nb3,nb3+1,nb3+2);
  max_volume_ = nb3[0]*nb3[1]*nb3[2];
  max_volume_ *= block_volume_[0];  // in James' code, this is block_volume_[max_level_ - min_level_]

  *volume_(block) = 0;

  
  compute_ (block);
}

//----------------------------------------------------------------------

double EnzoMethodMultipole::timestep ( Block * block ) throw()
{
  // initialize_(block);

  return timeStep_;
}

//======================================================================


void EnzoMethodMultipole::compute_ (Block * block) throw()
{
  EnzoBlock * enzo_block = enzo::block(block);
  int level = enzo_block->level();

  // double * mass = pmass(block);
  // double * com = pcom(block);
  // double * quadrupole = pquadrupole(block);
  
  if (block->is_leaf()) {
    
    compute_multipoles_ (block);
    begin_up_cycle_ (enzo_block);
    
  } else {

    if (0 <= level && level < cello::hierarchy()->max_level() ) {
      restrict_recv(enzo_block, nullptr);
    }
  }
    
}


void EnzoMethodMultipole::begin_up_cycle_(EnzoBlock * enzo_block) throw()
{
  const int level = enzo_block->level();

  double * mass = pmass(enzo_block);
  double * com = pcom(enzo_block);
  double * quadrupole = pquadrupole(enzo_block);

  if (level == 0) {

    CkPrintf("up cycle level: 0\n");
    CkPrintf("total mass: %f\n", *mass); 
    CkPrintf("COM: (%f, %f, %f)\n", com[0], com[1], com[2]);
    CkPrintf("quadrupole: ");
    for (int i = 0; i < 9; i++) {
      CkPrintf("%f  ", quadrupole[i]);
    }
    CkPrintf("\n");

  }

  else {

    // CkPrintf("up cycle level: %d\n", level);
    // CkPrintf("total mass: %f\n", *mass); 
    // CkPrintf("COM: (%f, %f, %f)\n", com[0], com[1], com[2]);
    // CkPrintf("quadrupole: ");
    // for (int i = 0; i < 9; i++) {
    //   CkPrintf("%f  ", quadrupole[i]);
    // }
    // CkPrintf("\n");

    restrict_send (enzo_block);

  }

  // enzo_block->compute_done();

  CkCallback callback(CkIndex_EnzoBlock::r_method_multipole_dualwalk_barrier(nullptr),
		      enzo::block_array());
  enzo::block(enzo_block)->contribute(callback);

}


void EnzoBlock::r_method_multipole_dualwalk_barrier(CkReductionMsg* msg)
{
  EnzoMethodMultipole * method =
    static_cast<EnzoMethodMultipole*> (this->method());

  if (cello::hierarchy()->max_level() == 0) {
    method->begin_down_cycle_(this);
  }

  else {
    method->dual_walk_(this);
  }
  

}

void EnzoMethodMultipole::dual_walk_(EnzoBlock * enzo_block) throw()
{

  Index index = enzo_block->index();
  int level = enzo_block->level();

  *volume_(enzo_block) = block_volume_[level];

  if (level == 0) {

    // call interact on all pairs of root blocks

    int type = enzo_block->is_leaf() ? 1 : 0;
    enzo::block_array()[index].p_method_multipole_traverse(index, type);
  }

  if (! enzo_block->is_leaf()) {
    
    // traverse terminates when all leaf blocks done. Since we use a
    // barrier on all blocks, we call the barrier now on all non-leaf
    // blocks
    CkCallback callback(CkIndex_EnzoBlock::r_method_multipole_traverse_complete(nullptr),
			enzo::block_array());
    enzo::block(enzo_block)->contribute(callback);
  }

}

void EnzoBlock::r_method_multipole_traverse_complete(CkReductionMsg * msg)
{
  // delete msg; is this necessary?

  EnzoMethodMultipole * method =
    static_cast<EnzoMethodMultipole*> (this->method());

  if (level() == 0){ //only the block at the top invokes begin_down_cycle
      method->begin_down_cycle_(this);
  }
}

void EnzoMethodMultipole::begin_down_cycle_(EnzoBlock * enzo_block) throw()
{
  const int level = enzo_block->level();

  double * c1 = pc1(enzo_block);
  double * c2 = pc2(enzo_block);
  double * c3 = pc3(enzo_block);

  if (enzo_block->is_leaf()) {

    CkPrintf("down cycle level: %d\n", level);

    CkPrintf("c1: ");
    for (int i = 0; i < 3; i++) {
      CkPrintf("%f  ", c1[i]);
    }
    CkPrintf("\n");

    CkPrintf("c2: ");
    for (int i = 0; i < 9; i++) {
      CkPrintf("%f  ", c2[i]);
    }
    CkPrintf("\n");

    CkPrintf("c3: ");
    for (int i = 0; i < 27; i++) {
      CkPrintf("%f  ", c3[i]);
    }
    CkPrintf("\n\n");

    evaluate_force_(enzo_block);

  }

  else {

    CkPrintf("down cycle level: %d\n", level);
    
    CkPrintf("c1: ");
    for (int i = 0; i < 3; i++) {
      CkPrintf("%f  ", c1[i]);
    }
    CkPrintf("\n");

    CkPrintf("c2: ");
    for (int i = 0; i < 9; i++) {
      CkPrintf("%f  ", c2[i]);
    }
    CkPrintf("\n\n");

    prolong_send (enzo_block);

  }

enzo_block->compute_done();
  
}


void EnzoMethodMultipole::restrict_send(EnzoBlock * enzo_block) throw()
{
  MultipoleMsg * msg = pack_multipole_(enzo_block);

  Index index_parent = enzo_block->index().index_parent(0);

  enzo::block_array()[index_parent].p_method_multipole_restrict_recv(msg);
}


void EnzoBlock::p_method_multipole_restrict_recv(MultipoleMsg * msg)
{
  EnzoMethodMultipole * method = 
    static_cast<EnzoMethodMultipole*> (this->method());

  method->restrict_recv(this,msg);
}


void EnzoMethodMultipole::restrict_recv
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{
  
  // Save multipole message from child
  if (msg != nullptr)
    {
      *pmsg_restrict(enzo_block, msg->child_index()) = msg;
    }
  
  // Continue if all expected messages received
  if ( psync_restrict(enzo_block)->next() ) {

    // Restore saved messages then clear
    for (int i = 0; i < cello::num_children(); i++) {
      msg = *pmsg_restrict(enzo_block, i);
      *pmsg_restrict(enzo_block, i) = nullptr;

      // Unpack multipoles from message then delete message
      unpack_multipole_(enzo_block, msg);
    }
  
    begin_up_cycle_ (enzo_block);
  }
  
}

void EnzoMethodMultipole::prolong_send(EnzoBlock * enzo_block) throw()
{

  ItChild it_child(cello::rank());
  int ic3[3];

  // std::string string = enzo_block->name();
  // CkPrintf("parent block name = %s\n", string.c_str());
    

  while (it_child.next(ic3)) {

    MultipoleMsg * msg = pack_coeffs_(enzo_block);

    Index index_child = enzo_block->index().index_child(ic3, cello::hierarchy()->min_level());

    //std::string string = enzo_block->name(index_child);
    //CkPrintf("sending to block name = %s\n", string.c_str());
    //fflush(stdout);

    enzo::block_array()[index_child].p_method_multipole_prolong_recv(msg);  

  }
}


void EnzoBlock::p_method_multipole_prolong_recv(MultipoleMsg * msg)
{
  EnzoMethodMultipole * method = 
    static_cast<EnzoMethodMultipole*> (this->method());

  method->prolong_recv(this, msg);
}


void EnzoMethodMultipole::prolong_recv
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{

  // Save multipole message from parent
  if (msg != nullptr)
    {

      MultipoleMsg ** temp = pmsg_prolong(enzo_block);

      if (temp == nullptr) {
        ERROR("early exit", "early exit");
      }
      *temp = msg;
    }

  //CkPrintf("message saved in prolong_recv\n");
  
  // Return if not ready yet 
  //if (! psync_prolong(enzo_block)->next() ) return;
  
  //CkPrintf("incremented sync counter in prolong_recv\n");

  // Restore saved message then clear
  msg = *pmsg_prolong(enzo_block);
  *pmsg_prolong(enzo_block) = nullptr;

  // Unpack coefficients from message then delete message
  unpack_coeffs_(enzo_block, msg);
    
  
  begin_down_cycle_ (enzo_block);
  
}


MultipoleMsg * EnzoMethodMultipole::pack_multipole_(EnzoBlock * enzo_block) throw()
{
 
  // Create a MultipoleMsg for sending data to parent
 
  MultipoleMsg * msg  = new MultipoleMsg;

  double * mass = pmass(enzo_block);
  double * com = pcom(enzo_block);
  double * quadrupole = pquadrupole(enzo_block);
 
  msg->mass = *mass;

  for (int i = 0; i < 3; i++) {
    msg->com[i] = com[i];   // use copy() ?
  }

  for (int i = 0; i < 9; i++) {
    msg->quadrupole[i] = quadrupole[i];
  }

  const int level = enzo_block->level();
  if (level > 0) {
    enzo_block->index().child(level,msg->ic3,msg->ic3+1,msg->ic3+2);
  }

  return msg;

}

// do I need int ic3[3] as an argument? 
MultipoleMsg * EnzoMethodMultipole::pack_coeffs_(EnzoBlock * enzo_block) throw()
{
 
  // Create a MultipoleMsg for sending data to children
 
  MultipoleMsg * msg  = new MultipoleMsg;

  double * com = pcom(enzo_block);
  double * c1 = pc1(enzo_block);
  double * c2 = pc2(enzo_block);
  double * c3 = pc3(enzo_block);

  for (int i = 0; i < 3; i++) {
    msg->com[i] = com[i];
  }

  for (int i = 0; i < 3; i++) {
    msg->c1[i] = c1[i];
  }

  for (int i = 0; i < 9; i++) {
    msg->c2[i] = c2[i];
  }
  
  for (int i = 0; i < 27; i++) {
    msg->c3[i] = c3[i];
  }

  return msg;

}


void EnzoMethodMultipole::unpack_multipole_
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{

  double * this_mass = pmass(enzo_block);
  double * this_com = pcom(enzo_block);
  double * this_quadrupole = pquadrupole(enzo_block);

  // copy data from msg to this EnzoBlock
  double child_mass = msg->mass;
  double * child_com = msg->com;
  double * child_quadrupole = msg->quadrupole;

  if (child_mass != 0) {

    double new_com[3];
    for (int i = 0; i < 3; i++) {
      new_com[i] = (*this_mass * this_com[i] + child_mass * child_com[i]) / (*this_mass + child_mass);
    }
    
    std::vector<double> shifted_this_quadrupole = shift_quadrupole_(this_quadrupole, *this_mass, this_com, new_com);
    std::vector<double> shifted_child_quadrupole = shift_quadrupole_(child_quadrupole, child_mass, child_com, new_com);
    for (int i = 0; i < 9; i++) {
        this_quadrupole[i] = shifted_this_quadrupole[i] + shifted_child_quadrupole[i];
    }

    for (int i = 0; i < 3; i++)
      this_com[i] = new_com[i];

    *this_mass += child_mass;
  }

  delete msg;

  // delete new_com ?
}

void EnzoMethodMultipole::unpack_coeffs_
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{
  // coefficient data from the current block
  std::vector<double> this_com (pcom(enzo_block), pcom(enzo_block) + 3);
  std::vector<double> this_c1  (pc1(enzo_block), pc1(enzo_block) + 3);
  std::vector<double> this_c2  (pc2(enzo_block), pc2(enzo_block) + 9);
  std::vector<double> this_c3  (pc3(enzo_block), pc3(enzo_block) + 27);
  
  // copy data from msg to this EnzoBlock
  std::vector<double> parent_com (msg->com, msg->com + 3);
  std::vector<double> parent_c1  (msg->c1, msg->c1 + 3);
  std::vector<double> parent_c2  (msg->c2, msg->c2 + 9);
  std::vector<double> parent_c3  (msg->c3, msg->c3 + 27);

  std::vector<double> com_shift = subtract_(parent_com, this_com, 3);

  std::vector<double> shifted_parent_c1_secondterm = dot_12_(com_shift, parent_c2);
  std::vector<double> shifted_parent_c1_thirdterm = dot_23_(outer_(com_shift, com_shift), dot_scalar_(0.5, parent_c3, 27));
  std::vector<double> shifted_parent_c1 = add_(add_(parent_c1, shifted_parent_c1_secondterm, 3), shifted_parent_c1_thirdterm, 3);
  
  std::vector<double> shifted_parent_c2_secondterm = dot_13_(com_shift, parent_c3);
  std::vector<double> shifted_parent_c2 = add_(parent_c2, shifted_parent_c2_secondterm, 9);
   
  std::vector<double> new_c1 = add_(this_c1, shifted_parent_c1, 3);   
  std::vector<double> new_c2 = add_(this_c2, shifted_parent_c2, 9);
  std::vector<double> new_c3 = add_(this_c3, parent_c3, 27);

  for (int i = 0; i < 3; i++) {
    pc1(enzo_block)[i] = new_c1[i];
  }

  for (int i = 0; i < 9; i++) {
    pc2(enzo_block)[i] = new_c2[i];
  }

  for (int i = 0; i < 27; i++) {
    pc3(enzo_block)[i] = new_c3[i];
  }

  delete msg;
}


/************************************************************************/
/// FMM functions

void EnzoMethodMultipole::compute_multipoles_ (Block * block) throw()
{
  Data * data = block->data();
  Field field = data->field();

  const int id_dens_ = field.field_id ("density");
  enzo_float * density = (enzo_float *) field.values("density");

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;

  field.dimensions  (id_dens_, &mx, &my, &mz);
  field.ghost_depth (id_dens_, &gx, &gy, &gz);

  block->cell_width(&hx, &hy, &hz);
  hz = 1.0; // this is to make math easier -- delete later
  double cell_vol = hx * hy * hz;
  // CkPrintf("Cell dims: (%f, %f, %f)\n", hx, hy, hz);

  double lo[3];
  block->lower(lo, lo+1, lo+2);

  double weighted_sum[3] = {};

  double * mass = pmass(block);
  double * com = pcom(block);
  double * quadrupole = pquadrupole(block);

  double cell_mass = 0;

  // should I have an if-statement with "rank" here to separate out 1D, 2D, and 3D?

  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {
        
        double dens = density[ix + mx * (iy + iz * my)];
        // CkPrintf("dens: %f\n", dens);

        cell_mass += dens;

        // CkPrintf("x coord: %f\n", (lo[0] + (ix-gx + 0.5)*hx));
        // CkPrintf("y coord: %f\n", (lo[1] + (iy-gy + 0.5)*hy));

	      weighted_sum[0] += dens * (lo[0] + (ix-gx + 0.5)*hx); 
	      weighted_sum[1] += dens * (lo[1] + (iy-gy + 0.5)*hy); 
	      weighted_sum[2] += dens * (lo[2] + (iz-gz + 0.5)*hz);

      }
    }
  }

  
  if (cell_mass != 0) {
    
    *mass += cell_mass; 

    com[0] = weighted_sum[0] / *mass;
    com[1] = weighted_sum[1] / *mass;
    com[2] = weighted_sum[2] / *mass;

    *mass *= cell_vol;

    for (int iz = gz; iz < mz-gz; iz++) {
      for (int iy = gy; iy < my-gy; iy++) {
        for (int ix = gx; ix < mx-gx; ix++) {
          
          double dens = density[ix + mx * (iy + iz * my)];

          double disp[3];
          disp[0] = (lo[0] + (ix-gx + 0.5)*hx) - com[0]; 
          disp[1] = (lo[1] + (iy-gy + 0.5)*hy) - com[1];
          disp[2] = (lo[2] + (iz-gz + 0.5)*hz) - com[2]; 

          for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
              quadrupole[3*i + j] += dens * cell_vol * disp[i] * disp[j];
            }
          }

        }
      }
    }
  }
}


void EnzoMethodMultipole::evaluate_force_(Block * block) throw()
{
  Data * data = block->data();
  Field field = data->field();
  
  enzo_float * dens = (enzo_float*) field.values ("density");
  enzo_float * accel_x = (enzo_float*) field.values ("acceleration_x");
  enzo_float * accel_y = (enzo_float*) field.values ("acceleration_y");
  enzo_float * accel_z = (enzo_float*) field.values ("acceleration_z");

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  block->cell_width(&hx, &hy, &hz);
  hz = 1.0; // this is to make math easier -- delete later
  double cell_vol = hx * hy * hz;

  double lo[3];
  block->lower(lo, lo+1, lo+2);

  std::vector<double> com (pcom(block), pcom(block) + 3);
  std::vector<double> c1 (pc1(block), pc1(block) + 3);
  std::vector<double> c2 (pc2(block), pc2(block) + 9);
  std::vector<double> c3 (pc3(block), pc3(block) + 27);
  
  // loop over all cells
  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {

        int i = ix + mx*(iy + my*iz);

        // compute force sourced from other Blocks
        std::vector<double> a (3, 0);
        a[0] = (lo[0] + (ix-gx + 0.5)*hx) - com[0]; 
        a[1] = (lo[1] + (iy-gy + 0.5)*hy) - com[1];
        a[2] = (lo[2] + (iz-gz + 0.5)*hz) - com[2];
        
        std::vector<double> second_term = dot_12_(a, c2);
        std::vector<double> third_term = dot_23_(outer_(a, a), dot_scalar_(0.5, c3, 27));

        // how does the code treat G?
        std::vector<double> block_force = add_(subtract_(c1, second_term, 3), third_term, 3);

        // subtracting rather than adding since block_force is multiplied by -G
        accel_x[i] -= block_force[0];
        accel_y[i] -= block_force[1];
        accel_z[i] -= block_force[2];

        CkPrintf("position: %f, %f, %f\n", (lo[0] + (ix-gx + 0.5)*hx), (lo[1] + (iy-gy + 0.5)*hy), (lo[2] + (iz-gz + 0.5)*hz));
        CkPrintf("Block force: %f, %f, %f\n", -1.0*block_force[0], -1.0*block_force[1], -1.0*block_force[2]);

        std::vector<double> tot_cell_force (3, 0); // for debugging purposes

        // compute force sourced from other cells
        for (int iz2 = gz; iz2 < mz-gz; iz2++) {
          for (int iy2 = gy; iy2 < my-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx-gx; ix2++) {

              int i2 = ix2 + mx*(iy2 + my*iz2);

              if (i != i2) {

                // displacement vector pointing from current cell to interacting cell
                std::vector<double> disp (3, 0);
                disp[0] = (ix2 - ix) * hx;
                disp[1] = (iy2 - iy) * hy;
                disp[2] = (iz2 - iz) * hz;

                //CkPrintf("Cell mass: %f\n", dens[i2]*cell_vol);
                std::vector<double> cell_force = newton_force_(dens[i2]*cell_vol, disp); 

                accel_x[i] += cell_force[0];
                accel_y[i] += cell_force[1];
                accel_z[i] += cell_force[2];

                // CkPrintf("Cell 2: %d\n", i2);
                // CkPrintf("disp: %f, %f, %f\n", disp[0], disp[1], disp[2]);
                // CkPrintf("Cell force: %f, %f, %f\n", cell_force[0], cell_force[1], cell_force[2]);

                tot_cell_force[0] += cell_force[0];
                tot_cell_force[1] += cell_force[1];
                tot_cell_force[2] += cell_force[2];
              }
            }
          }
        }

        CkPrintf("Cell force: %f, %f, %f\n", tot_cell_force[0], tot_cell_force[1], tot_cell_force[2]);
        CkPrintf("Accel: %f, %f, %f\n\n", accel_x[i], accel_y[i], accel_z[i]);
      }
    }
  }
} 

/************************************************************************/

/// Dual Tree Walk helper functions

void EnzoBlock::p_method_multipole_traverse(Index index, int type)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->traverse(this, index, type);
}


void EnzoMethodMultipole::traverse
(EnzoBlock * enzo_block, Index index_b, int type_b)
{
  Index index_a = enzo_block->index();
  const bool known_leaf_b = (type_b != -1);
  const bool is_leaf_a = enzo_block->is_leaf();

  if (! known_leaf_b) {
    // If unknown whether B is leaf or not, flip arguments and send to B
    enzo::block_array()[index_b].p_method_multipole_traverse
      (index_a, is_leaf_a ? 1 : 0);
    return;
  }

  // Determine if blocks are "far", and save radii for later
  double ra,rb;
  bool mac = is_far_(enzo_block,index_b,&ra,&rb);

  // is_leaf values: 0 no, 1 yes -1 unknown
  const bool is_leaf_b = (type_b == 1);

  ItChild it_child_a(cello::rank());
  ItChild it_child_b(cello::rank());
  int ica3[3],icb3[3];

  int level_a = index_a.level();
  int level_b = index_b.level();
  int volume_a = block_volume_[level_a];
  int volume_b = block_volume_[level_b];

  if (index_a == index_b) {

    // loop 1 through A children
    while (it_child_a.next(ica3)) {
      int ica = ica3[0] + 2*(ica3[1] + 2*ica3[2]);
      Index index_a_child = index_a.index_child(ica3, 0);

      // loop 2 through A children
      while (it_child_b.next(icb3)) {
        int icb = icb3[0] + 2*(icb3[1] + 2*icb3[2]);
        Index index_b_child = index_b.index_child(icb3, 0);
        // avoid calling both traverse (A,B) and traverse (B,A)

        if (ica <= icb) {  
          // call traverse on child 1 and child 2
          enzo::block_array()[index_a_child].p_method_multipole_traverse
            (index_b_child,-1);
        }
      }
    }
  }

  else if (is_leaf_a && is_leaf_b) {

    // interact
    traverse_direct_pair (enzo_block,index_a,volume_a,index_b,volume_b);

  }

  else if (mac) {

    // interact
    traverse_approx_pair (enzo_block,index_a,volume_a,index_b,volume_b);

  } 
  
  
  else if (is_leaf_a) {

    while (it_child_b.next(icb3)) {

      Index index_b_child = index_b.index_child(icb3, 0);
      traverse(enzo_block,index_b_child,-1);
    }

  } 
  
  else if (is_leaf_b) {

    while (it_child_a.next(ica3)) {

      Index index_a_child = index_a.index_child(ica3, 0);
      enzo::block_array()[index_a_child].p_method_multipole_traverse
        (index_b,type_b);
    }
  } 

  else if (ra < rb) { // open the larger block (b)

    while (it_child_b.next(icb3)) {
      Index index_b_child = index_b.index_child(icb3, 0);
      traverse(enzo_block,index_b_child,-1);
    }
  } 
  
  else {  // open the larger block (a)

    while (it_child_a.next(ica3)) {
      Index index_a_child = index_a.index_child(ica3, 0);
      enzo::block_array()[index_a_child].p_method_multipole_traverse
        (index_b,type_b);   
    }
  }
}


void EnzoMethodMultipole::traverse_approx_pair
(EnzoBlock * enzo_block,
 Index index_a, int volume_a,
 Index index_b, int volume_b)
{

  // compute the a->b interaction
  MultipoleMsg * msg_a = pack_multipole_(enzo_block);
  enzo::block_array()[index_b].p_method_multipole_interact_approx(msg_a);

  // compute the b->a interaction
  enzo::block_array()[index_b].p_method_multipole_interact_approx_send (index_a);

  // update volumes for each block to detect when to terminate traverse
  
  update_volume (enzo_block,index_b,volume_b);
  // do we need to check if blocks are equal?
  if (index_a != index_b) {
    enzo::block_array()[index_b].p_method_multipole_update_volume (index_a,volume_a);
  }
}

void EnzoBlock::p_method_multipole_interact_approx (MultipoleMsg * msg)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->interact_approx_ (this, msg);
}

void EnzoBlock::p_method_multipole_interact_approx_send (Index receiver)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->interact_approx_send (this, receiver);
}

void EnzoMethodMultipole::interact_approx_send(EnzoBlock * enzo_block, Index receiver) throw()
{
  MultipoleMsg * msg = pack_multipole_(enzo_block);
  enzo::block_array()[receiver].p_method_multipole_interact_approx(msg);
}

void EnzoMethodMultipole::interact_approx_(Block * block, MultipoleMsg * msg_b) throw()
{
  
  std::vector<double> com_a (pcom(block), pcom(block) + 3);
  std::vector<double> c1_a  (pc1(block), pc1(block) + 3);
  std::vector<double> c2_a  (pc2(block), pc2(block) + 9);
  std::vector<double> c3_a  (pc3(block), pc3(block) + 27);

  double mass_b = msg_b->mass;
  std::vector<double> com_b (msg_b->com, msg_b->com + 3);
  std::vector<double> quadrupole_b (msg_b->quadrupole, msg_b->quadrupole + 9);


  std::vector<double> rvec = subtract_(com_b, com_a, 3);         // displacement vector between com_b and com_a
  double r = sqrt(dot_11_(rvec, rvec));                          // magnitude of displacement vector

  std::vector<double> d1 = dot_scalar_(-1.0/pow(r,3), rvec, 3);  // derivative tensor d1
  std::vector<double> d2 (9, 0);                                 // derivative tensor d2
  std::vector<double> d3 (27, 0);                                // derivative tensor d3
       
          
  // compute the components of d2
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {

      d2[3*i + j] = 3.0/pow(r,5) * rvec[i] * rvec[j];

      if (i == j) {
        d2[3*i + j] -= 1.0/pow(r,3);
      }
    }
  }
      
  // compute the components of d3
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {

        d3[9*k + 3*i + j] = -15.0/pow(r,7) * rvec[i] * rvec[j] * rvec[k];

        if (i == j) {
          d3[9*k + 3*i + j] += 3.0/pow(r,5) * rvec[k];
        }

        if (j == k) {
          d3[9*k + 3*i + j] += 3.0/pow(r,5) * rvec[i];
        }

        if (i == k) {
          d3[9*k + 3*i + j] += 3.0/pow(r,5) * rvec[j];
        }
      }
    }
  }
     
     
  // compute the coefficients of the Taylor expansion of acceleration due to the particles in Block b
  std::vector<double> delta_c1 = add_(dot_scalar_(mass_b, d1, 3), dot_23_(quadrupole_b, dot_scalar_(0.5, d3, 27)), 3);
  std::vector<double> delta_c2 = dot_scalar_(mass_b, d2, 9);
  std::vector<double> delta_c3 = dot_scalar_(mass_b, d3, 27);
   
  // add the coefficients for the new interaction to the coefficients already associated with this Block
  std::vector<double> new_c1 = add_(c1_a, delta_c1, 3);  
  std::vector<double> new_c2 = add_(c2_a, delta_c2, 9);
  std::vector<double> new_c3 = add_(c3_a, delta_c3, 27);

  for (int i = 0; i < 3; i++) {
    pc1(block)[i] = new_c1[i];
  }

  for (int i = 0; i < 9; i++) {
    pc2(block)[i] = new_c2[i];
  }

  for (int i = 0; i < 27; i++) {
    pc3(block)[i] = new_c3[i];
  }

}



void EnzoMethodMultipole::traverse_direct_pair
(EnzoBlock * enzo_block,
 Index index_a, int volume_a,
 Index index_b, int volume_b)
{

  // compute the a->b interaction
  pack_dens_(enzo_block, index_b);

  // compute the b->a interaction
  enzo::block_array()[index_b].p_method_multipole_interact_direct_send (index_a);

  // update volumes for each block to detect when to terminate traverse
  update_volume (enzo_block,index_b,volume_b);
  if (index_a != index_b) {
    // (note blocks may be equal, so don't double-count)
    enzo::block_array()[index_b].p_method_multipole_update_volume (index_a,volume_a);
  }

}

void EnzoBlock::p_method_multipole_interact_direct (int n, char * msg)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->interact_direct_ (this, msg);
}

void EnzoBlock::p_method_multipole_interact_direct_send (Index receiver)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->interact_direct_send (this, receiver);
}

void EnzoMethodMultipole::interact_direct_send(EnzoBlock * enzo_block, Index receiver) throw()
{
  pack_dens_(enzo_block, receiver);
}


void EnzoMethodMultipole::pack_dens_(EnzoBlock * enzo_block, Index index_b) throw()
{
  Data * data = enzo_block->data();
  Field field = data->field();

  enzo_float * dens = (enzo_float*) field.values ("density");

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;
  double lo[3];

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  enzo_block->cell_width(&hx, &hy, &hz);
  hz = 1.0; // this is to make math easier -- delete later

  enzo_block->lower(lo, lo+1, lo+2);

  int size = 0;

  SIZE_SCALAR_TYPE(size, int, mx);
  SIZE_SCALAR_TYPE(size, int, my);
  SIZE_SCALAR_TYPE(size, int, mz);
  SIZE_SCALAR_TYPE(size, double, hx);
  SIZE_SCALAR_TYPE(size, double, hy);
  SIZE_SCALAR_TYPE(size, double, hz);
  SIZE_ARRAY_TYPE(size, double, lo, 3);
  SIZE_ARRAY_TYPE(size, enzo_float, dens, mx*my*mz);

  char * buffer = new char[size];

  char * pc; 
  pc = buffer;

  SAVE_SCALAR_TYPE(pc, int, mx);
  SAVE_SCALAR_TYPE(pc, int, my);
  SAVE_SCALAR_TYPE(pc, int, mz);
  SAVE_SCALAR_TYPE(pc, double, hx);
  SAVE_SCALAR_TYPE(pc, double, hy);
  SAVE_SCALAR_TYPE(pc, double, hz);
  SAVE_ARRAY_TYPE(pc, double, lo, 3);
  SAVE_ARRAY_TYPE(pc, enzo_float, dens, mx*my*mz);


  enzo::block_array()[index_b].p_method_multipole_interact_direct (size, buffer);
}


void EnzoMethodMultipole::interact_direct_(Block * block, char * msg_b) throw()
{

  Data * data = block->data();
  Field field = data->field();

  enzo_float * accel_x = (enzo_float*) field.values ("acceleration_x");
  enzo_float * accel_y = (enzo_float*) field.values ("acceleration_y");
  enzo_float * accel_z = (enzo_float*) field.values ("acceleration_z");


  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  block->cell_width(&hx, &hy, &hz);
  hz = 1.0; // this is to make math easier -- delete later
  double cell_vol = hx * hy * hz;

  double lo[3];
  block->lower(lo, lo+1, lo+2);

  // unpacking data from msg_b
  int mx2, my2, mz2;
  double hx2, hy2, hz2;
  double lo2[3];

  char * pc;
  pc = msg_b;

  LOAD_SCALAR_TYPE(pc, int, mx2);
  LOAD_SCALAR_TYPE(pc, int, my2);
  LOAD_SCALAR_TYPE(pc, int, mz2);
  LOAD_SCALAR_TYPE(pc, double, hx2);
  LOAD_SCALAR_TYPE(pc, double, hy2);
  LOAD_SCALAR_TYPE(pc, double, hz2);
  LOAD_ARRAY_TYPE(pc, double, lo2, 3);

  double dens[mx2*my2*mz2];
  LOAD_ARRAY_TYPE(pc, double, dens, mx2*my2*mz2);


  // loop over all cells in this Block
  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {

        int i = ix + mx*(iy + my*iz);

        CkPrintf("position: %f, %f, %f\n", (lo[0] + (ix-gx + 0.5)*hx), (lo[1] + (iy-gy + 0.5)*hy), (lo[2] + (iz-gz + 0.5)*hz));

        // compute force sourced from cells in Block b
        for (int iz2 = gz; iz2 < mz2-gz; iz2++) {
          for (int iy2 = gy; iy2 < my2-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx2-gx; ix2++) {

              int i2 = ix2 + mx2*(iy2 + my2*iz2);

              // disp points from cell in current Block to cell in Block b
              std::vector<double> disp (3, 0);
              disp[0] = lo2[0] - lo[0] + (ix2 - ix)*hx; 
              disp[1] = lo2[1] - lo[1] + (iy2 - iy)*hy;
              disp[2] = lo2[2] - lo[2] + (iz2 - iz)*hz;
                
              std::vector<double> cell_force = newton_force_(dens[i2]*cell_vol, disp); 

              accel_x[i] += cell_force[0];
              accel_y[i] += cell_force[1];
              accel_z[i] += cell_force[2];

            }
          }
        }

        // CkPrintf("Leaf-leaf cell force: %f, %f, %f\n\n", accel_x[i], accel_y[i], accel_z[i]);
      }
    }
  }
}


void EnzoBlock::p_method_multipole_update_volume (Index index, int volume)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->update_volume (this, index, volume);
}

// i don't think we need the index?
void EnzoMethodMultipole::update_volume
(EnzoBlock * enzo_block, Index index, int volume)
{

  std::string string = enzo_block->name();

  if (enzo_block->is_leaf()) {

    *volume_(enzo_block) += volume;

    if (*volume_(enzo_block) == max_volume_) {

      CkCallback callback(CkIndex_EnzoBlock::r_method_multipole_traverse_complete(nullptr),
			        enzo::block_array());
      enzo_block->contribute(callback);

    }

  } 
  
  else {

    const int min_level = 0; //cello::hierarchy()->min_level();
    ItChild it_child(cello::rank());
    int ic3[3];

    while (it_child.next(ic3)) {

      Index index_child = enzo_block->index().index_child(ic3,min_level);
      enzo::block_array()[index_child].p_method_multipole_update_volume (index,volume);
    }
  }
}

bool EnzoMethodMultipole::is_far_ (EnzoBlock * enzo_block,
                             Index index_b, double *ra, double *rb) const
{
  double iam3[3],iap3[3],ibm3[3],ibp3[3],na3[3],nb3[3];
  enzo_block->lower(iam3,iam3+1,iam3+2);
  enzo_block->upper(iap3,iap3+1,iap3+2);
  enzo_block->lower(ibm3,ibm3+1,ibm3+2,&index_b);  // overloaded lower/upper? (see Cello/mesh_Block.hpp)
  enzo_block->upper(ibp3,ibp3+1,ibp3+2,&index_b);
  *ra=0;
  *rb=0;
  double ra3[3]={0,0,0},rb3[3]={0,0,0},d3[3]={0,0,0};
  double d = 0;

  for (int i=0; i<cello::rank(); i++) {
    ra3[i] = (iap3[i] - iam3[i]);
    rb3[i] = (ibp3[i] - ibm3[i]);
    *ra += ra3[i]*ra3[i];
    *rb += rb3[i]*rb3[i];
    double ca = (iam3[i]+0.5*ra3[i]);
    double cb = (ibm3[i]+0.5*rb3[i]);
    d3[i] = (ca - cb);
    d += d3[i]*d3[i];
  }

  d = sqrt(d);
  *ra = 0.5*sqrt(*ra);
  *rb = 0.5*sqrt(*rb);

  return  (d * theta_ > (*ra + *rb));
}

/************************************************************************/