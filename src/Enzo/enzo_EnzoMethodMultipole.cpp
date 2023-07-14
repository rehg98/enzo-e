// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMultipole.cpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     July 19, 2022
/// @brief    Compute multipoles and pass multipoles up octree

#include "cello.hpp"

#include "enzo.hpp"

// #include <fstream> // necessary?

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

  // is there a cleaner way of zeroing out everything?
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
  

  Data * data = block->data();
  Field field = data->field();

  enzo_float * dens = (enzo_float *) field.values("density");
  enzo_float * accel_x = (enzo_float *) field.values("acceleration_x");
  enzo_float * accel_y = (enzo_float *) field.values("acceleration_y");
  enzo_float * accel_z = (enzo_float *) field.values("acceleration_z");

  // CkPrintf("density: %p\n", dens);
  // CkPrintf("acceleration: %p\n", accel_x);

  // int mx, my, mz;
  // int gx, gy, gz;

  // field.dimensions  (0, &mx, &my, &mz);
  // field.ghost_depth (0, &gx, &gy, &gz);

  // std::string string = block->name();
  // CkPrintf("%s in compute\n", string.c_str());
  // for (int iz = gz; iz < mz-gz; iz++) {
  //   for (int iy = gy; iy < my-gy; iy++) {
  //     for (int ix = gx; ix < mx-gx; ix++) {
        
  //       int i = ix + mx * (iy + iz * my);
  //       CkPrintf("%f ", accel_x[i]);
  //       // accel_x[i] = 0;
  //       // accel_y[i] = 0;
  //       // accel_z[i] = 0;

  //     }
  //   }
  // }
  // CkPrintf("\n\n");




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


  if (block->is_leaf()) {

    // CkPrintf("a leaf in compute\n");

    int nprtls = 4;
    // mass, x, y, z, ax, ay, az
    double prtls[nprtls][7] = {{1.0, 0.3, 0.3, 0.5, 0, 0, 0},
                               {1.0, -0.3, -0.42, 0.5, 0, 0, 0},
                               {0.5, 0.35, 0.45, 0.5, 0, 0, 0},
                               {10.0, -0.5, 0.0, 0.5, 0, 0, 0}};
    InitializeParticles(block, nprtls, prtls);
  }
  
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

  // method->begin_down_cycle_(this);

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

    // std::string string = enzo_block->name(index_child);
    // CkPrintf("sending to block name = %s\n", string.c_str());
    // fflush(stdout);

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
    
    std::array<double, 9> shifted_this_quadrupole = shift_quadrupole_(this_quadrupole, *this_mass, this_com, new_com);
    std::array<double, 9> shifted_child_quadrupole = shift_quadrupole_(child_quadrupole, child_mass, child_com, new_com);
    for (int i = 0; i < 9; i++) {
        this_quadrupole[i] = shifted_this_quadrupole[i] + shifted_child_quadrupole[i];
    }

    for (int i = 0; i < 3; i++)
      this_com[i] = new_com[i];

    *this_mass += child_mass;

    // delete new_com ?
  }

  delete msg;
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

  std::string string = block->name();
  // const int id_dens_ = field.field_id ("density");
  enzo_float * density = (enzo_float *) field.values("density");
  // enzo_float * accel_x = (enzo_float *) field.values("acceleration_x");
  //CkPrintf("dens[2], block in compute_multipole: %f, %s\n", density[2], string.c_str());

  

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  block->cell_width(&hx, &hy, &hz);
  //hz = 1.0; // this is to make math easier -- delete later
  double cell_vol = hx * hy * hz;
  // CkPrintf("Cell dims: (%f, %f, %f)\n", hx, hy, hz);

  double lo[3];
  block->lower(lo, lo+1, lo+2);

  // CkPrintf("%s\n", string.c_str());
  // for (int i = 0; i < mx*my*mz; i++) {
  //   CkPrintf("%f ", density[i]);
  // }
  
  //std::string string = block->name();
  // CkPrintf("%s accel in compute_multipole\n", string.c_str());
  // for (int iz = gz; iz < mz-gz; iz++) {
  //   for (int iy = gy; iy < my-gy; iy++) {
  //     for (int ix = gx; ix < mx-gx; ix++) {
        
  //       double acc = accel_x[ix + mx * (iy + iz * my)];

  //       CkPrintf("(%d, %d, %d): %f  ", ix, iy, iz, acc);
  //     }
  //   }
  // }
  // CkPrintf("\n\n");


  Particle particle = block->data()->particle();
  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();
  const int num_is_grav = particle_groups->size("is_gravitating");

  // distinguish between mass as attribute vs. mass as constant??
  enzo_float * prtmass = NULL;
  int dm;



  double weighted_sum[3] = {};

  double * mass = pmass(block);
  double * com = pcom(block);
  double * quadrupole = pquadrupole(block);


  // can compute quadrupole in the same loop as the mass and com -- just need running sum over positions squared
  // ()
  // should I have an if-statement with "rank" here to separate out 1D, 2D, and 3D?

  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {
        
        double dens = density[ix + mx * (iy + iz * my)];

        // CkPrintf("(%d, %d, %d): %f  ", ix, iy, iz, dens);

        *mass += dens * cell_vol;

	      weighted_sum[0] += dens * cell_vol * (lo[0] + (ix-gx + 0.5)*hx); 
	      weighted_sum[1] += dens * cell_vol * (lo[1] + (iy-gy + 0.5)*hy); 
	      weighted_sum[2] += dens * cell_vol * (lo[2] + (iz-gz + 0.5)*hz);

      }
    }
  }

  

  // loop over particles, adding masses to *mass and adding mass*position to weighted_sum
  for (int ipt = 0; ipt < num_is_grav; ipt++) {
    const int it = particle.type_index(particle_groups->item("is_gravitating",ipt));

    int imass = 0;

    // check correct precision for position
    int ia = particle.attribute_index(it,"x");
    int ba = particle.attribute_bytes(it,ia); // "bytes (actual)"
    int be = sizeof(enzo_float);                // "bytes (expected)"

      ASSERT4 ("EnzoMethodPmUpdate::compute()",
	       "Particle type %s attribute %s defined as %s but expecting %s",
	       particle.type_name(it).c_str(),
	       particle.attribute_name(it,ia).c_str(),
	       ((ba == 4) ? "single" : ((ba == 8) ? "double" : "quadruple")),
	       ((be == 4) ? "single" : ((be == 8) ? "double" : "quadruple")),
	       (ba == be));

    for (int ib = 0; ib < particle.num_batches(it); ib++) {

      const int np = particle.num_particles(it,ib);

      if (particle.has_attribute(it,"mass")) {

        // Particle type has an attribute called "mass".
        // In this case we set prtmass to point to the mass attribute array
        // Also set dm to be the stride for the "mass" attribute
        imass = particle.attribute_index(it,"mass");
        prtmass = (enzo_float *) particle.attribute_array( it, imass, ib);
        dm = particle.stride(it,imass);

	    } 
      else {

        // Particle type has a constant called "mass".
        // In this case we set prtmass to point to the value
        // of the mass constant.
        // dm is set to 0, which will mean that we can loop through an
        // "array" of length 1.
        imass = particle.constant_index(it,"mass");
        prtmass = (enzo_float*)particle.constant_value(it,imass);
        dm = 0;
      }

      const int ia_x  = particle.attribute_index(it,"x");
      const int ia_y  = particle.attribute_index(it,"y");
      const int ia_z  = particle.attribute_index(it,"z");

      enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
      enzo_float * ya =  (enzo_float *)particle.attribute_array (it,ia_y,ib);
      enzo_float * za =  (enzo_float *)particle.attribute_array (it,ia_z,ib);

      const int dx =  particle.stride(it,ia_x);
      const int dy =  particle.stride(it,ia_y);
      const int dz =  particle.stride(it,ia_z);

      for (int ip=0; ip < np; ip++) {
        
        *mass += prtmass[ip*dm];

        // how are particle x, y, z defined? where is origin?
        weighted_sum[0] += prtmass[ip*dm] * xa[ip*dx]; 
	      weighted_sum[1] += prtmass[ip*dm] * ya[ip*dy]; 
	      weighted_sum[2] += prtmass[ip*dm] * za[ip*dz];

      }
    }
  }
  

  if (*mass != 0) { 

    com[0] = weighted_sum[0] / *mass;
    com[1] = weighted_sum[1] / *mass;
    com[2] = weighted_sum[2] / *mass;

    // can compute quadrupole in one loop?
    for (int iz = gz; iz < mz-gz; iz++) {
      for (int iy = gy; iy < my-gy; iy++) {
        for (int ix = gx; ix < mx-gx; ix++) {
          
          double dens = density[ix + mx * (iy + iz * my)];

          double disp[3]; // do I need to delete disp later?
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


    for (int ipt = 0; ipt < num_is_grav; ipt++) {
      const int it = particle.type_index(particle_groups->item("is_gravitating",ipt));

      int imass = 0;

      for (int ib = 0; ib < particle.num_batches(it); ib++) {

        const int np = particle.num_particles(it,ib);

        if (particle.has_attribute(it,"mass")) {
          imass = particle.attribute_index(it,"mass");
          prtmass = (enzo_float *) particle.attribute_array( it, imass, ib);
          dm = particle.stride(it,imass);
        } 
        else {
          imass = particle.constant_index(it,"mass");
          prtmass = (enzo_float*)particle.constant_value(it,imass);
          dm = 0;
        }

        const int ia_x  = particle.attribute_index(it,"x");
        const int ia_y  = particle.attribute_index(it,"y");
        const int ia_z  = particle.attribute_index(it,"z");

        enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
        enzo_float * ya =  (enzo_float *)particle.attribute_array (it,ia_y,ib);
        enzo_float * za =  (enzo_float *)particle.attribute_array (it,ia_z,ib);

        const int dx =  particle.stride(it,ia_x);
        const int dy =  particle.stride(it,ia_y);
        const int dz =  particle.stride(it,ia_z);

        for (int ip=0; ip < np; ip++) {

          double disp[3]; // do I need to delete disp later?
          disp[0] = xa[ip*dx] - com[0]; 
          disp[1] = ya[ip*dy] - com[1];
          disp[2] = za[ip*dz] - com[2]; 

          for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
              quadrupole[3*i + j] += prtmass[ip*dm] * disp[i] * disp[j];
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
  
  
  enzo_float * accel_x = (enzo_float*) field.values ("acceleration_x");
  enzo_float * accel_y = (enzo_float*) field.values ("acceleration_y");
  enzo_float * accel_z = (enzo_float*) field.values ("acceleration_z");

  enzo_float * dens = (enzo_float*) field.values ("density");

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  // std::string string = block->name();
  // CkPrintf("%s accel\n", string.c_str());
  // for (int iz = gz; iz < mz-gz; iz++) {
  //   for (int iy = gy; iy < my-gy; iy++) {
  //     for (int ix = gx; ix < mx-gx; ix++) {
        
  //       double acc = accel_x[ix + mx * (iy + iz * my)];

  //       CkPrintf("(%d, %d, %d): %f  ", ix, iy, iz, acc);
  //     }
  //   }
  // }
  // CkPrintf("\n\n");

  block->cell_width(&hx, &hy, &hz);
  // hz = 1.0; // this is to make math easier -- delete later
  double cell_vol = hx * hy * hz;

  double lo[3];
  block->lower(lo, lo+1, lo+2);


  Particle particle = block->data()->particle();
  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();
  const int num_is_grav = particle_groups->size("is_gravitating");
  const int num_prtl_types = particle_descr->num_types();

  // distinguish between mass as attribute vs. mass as constant??
  enzo_float * prtmass2 = NULL;
  int dm2;

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
                // consider replacing with array? (will this work with newton_force?)
                std::vector<double> disp (3, 0);
                disp[0] = (ix2 - ix) * hx;
                disp[1] = (iy2 - iy) * hy;
                disp[2] = (iz2 - iz) * hz;

                // CkPrintf("Cell dens: %f\n", dens[i2]);
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

                // CkPrintf("Partial cell force: %f, %f, %f\n", cell_force[0], cell_force[1], cell_force[2]);
              }
            }
          }
        }
        CkPrintf("Cell force: %f, %f, %f\n", tot_cell_force[0], tot_cell_force[1], tot_cell_force[2]);

        std::vector<double> tot_prt_force (3, 0); // for debugging purposes

        // compute cell force sourced from particles
        for (int ipt = 0; ipt < num_is_grav; ipt++) {
          const int it = particle.type_index(particle_groups->item("is_gravitating",ipt));

          int imass2 = 0;

          for (int ib = 0; ib < particle.num_batches(it); ib++) {

            const int np = particle.num_particles(it,ib);

            if (particle.has_attribute(it,"mass")) {
              imass2 = particle.attribute_index(it,"mass");
              prtmass2 = (enzo_float *) particle.attribute_array( it, imass2, ib);
              dm2 = particle.stride(it,imass2);
            } 
            else {
              imass2 = particle.constant_index(it,"mass");
              prtmass2 = (enzo_float*)particle.constant_value(it,imass2);
              dm2 = 0;
            }

            const int ia_x2  = particle.attribute_index(it,"x");
            const int ia_y2  = particle.attribute_index(it,"y");
            const int ia_z2  = particle.attribute_index(it,"z");

            enzo_float * xa2 =  (enzo_float *)particle.attribute_array (it,ia_x2,ib);
            enzo_float * ya2 =  (enzo_float *)particle.attribute_array (it,ia_y2,ib);
            enzo_float * za2 =  (enzo_float *)particle.attribute_array (it,ia_z2,ib);

            const int dx2 =  particle.stride(it,ia_x2);
            const int dy2 =  particle.stride(it,ia_y2);
            const int dz2 =  particle.stride(it,ia_z2);

            for (int ip=0; ip < np; ip++) {
              
              // disp points from current cell to particle
              std::vector<double> disp (3, 0);
              disp[0] = xa2[ip*dx2] - (lo[0] + (ix-gx + 0.5)*hx); 
              disp[1] = ya2[ip*dy2] - (lo[1] + (iy-gy + 0.5)*hy);
              disp[2] = za2[ip*dz2] - (lo[2] + (iz-gz + 0.5)*hz);
                
              std::vector<double> prtcell_force = newton_force_(prtmass2[ip*dm2], disp); 

              accel_x[i] += prtcell_force[0];
              accel_y[i] += prtcell_force[1];
              accel_z[i] += prtcell_force[2];

              tot_prt_force[0] += prtcell_force[0];
              tot_prt_force[1] += prtcell_force[1];
              tot_prt_force[2] += prtcell_force[2];
            }
          }
        }

        CkPrintf("Particle force: %f, %f, %f\n", tot_prt_force[0], tot_prt_force[1], tot_prt_force[2]);
        CkPrintf("Accel: %f, %f, %f\n\n", accel_x[i], accel_y[i], accel_z[i]);
      }
    }
  }


  //std::string string = block->name();
  // CkPrintf("%s dens\n", string.c_str());
  // for (int iz = gz; iz < mz-gz; iz++) {
  //   for (int iy = gy; iy < my-gy; iy++) {
  //     for (int ix = gx; ix < mx-gx; ix++) {
        
  //       double denses = dens[ix + mx * (iy + iz * my)];

  //       CkPrintf("(%d, %d, %d): %f  ", ix, iy, iz, denses);
  //     }
  //   }
  // }
  // CkPrintf("\n\n");
  // 



  // loop over all particles
  // change this to loop over all particles, not just the ones that are gravitating

  for (int it = 0; it < num_prtl_types; it++) {
    
    for (int ib = 0; ib < particle.num_batches(it); ib++) {

      const int np = particle.num_particles(it,ib);

      const int ia_x  = particle.attribute_index(it,"x");
      const int ia_y  = particle.attribute_index(it,"y");
      const int ia_z  = particle.attribute_index(it,"z");
      const int ia_ax  = particle.attribute_index(it,"ax");
      const int ia_ay  = particle.attribute_index(it,"ay");
      const int ia_az  = particle.attribute_index(it,"az");

      enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
      enzo_float * ya =  (enzo_float *)particle.attribute_array (it,ia_y,ib);
      enzo_float * za =  (enzo_float *)particle.attribute_array (it,ia_z,ib);
      enzo_float * axa =  (enzo_float *)particle.attribute_array (it,ia_ax,ib);
      enzo_float * aya =  (enzo_float *)particle.attribute_array (it,ia_ay,ib);
      enzo_float * aza =  (enzo_float *)particle.attribute_array (it,ia_az,ib);

      const int dx =  particle.stride(it,ia_x);
      const int dy =  particle.stride(it,ia_y);
      const int dz =  particle.stride(it,ia_z);
      const int dax =  particle.stride(it,ia_ax);
      const int day =  particle.stride(it,ia_ay);
      const int daz =  particle.stride(it,ia_az);

      for (int ip=0; ip < np; ip++) {

        // compute force sourced by other blocks
        std::vector<double> a (3, 0);
        a[0] =  xa[ip*dx] - com[0]; 
        a[1] =  ya[ip*dy] - com[1];
        a[2] =  za[ip*dz] - com[2];
        
        std::vector<double> second_term = dot_12_(a, c2);
        std::vector<double> third_term = dot_23_(outer_(a, a), dot_scalar_(0.5, c3, 27));

        // how does the code treat G?
        std::vector<double> block_force = add_(subtract_(c1, second_term, 3), third_term, 3);

        // subtracting rather than adding since block_force is multiplied by -G
        axa[ip*dax] -= block_force[0];
        aya[ip*day] -= block_force[1];
        aza[ip*daz] -= block_force[2];

        CkPrintf("position (prt): %f, %f, %f\n", xa[ip*dx], ya[ip*dy], za[ip*dz]);
        CkPrintf("Block force (prt): %f, %f, %f\n", -1.0*block_force[0], -1.0*block_force[1], -1.0*block_force[2]);

        std::vector<double> tot_cell_force (3, 0); // for debugging purposes
        
        // compute force sourced by cells
        for (int iz2 = gz; iz2 < mz-gz; iz2++) {
          for (int iy2 = gy; iy2 < my-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx-gx; ix2++) {

              int i2 = ix2 + mx*(iy2 + my*iz2);

              // disp points from current particle to interacting cell
              std::vector<double> disp (3, 0);
              disp[0] = (lo[0] + (ix2-gx + 0.5)*hx) - xa[ip*dx]; 
              disp[1] = (lo[1] + (iy2-gy + 0.5)*hy) - ya[ip*dy];
              disp[2] = (lo[2] + (iz2-gz + 0.5)*hz) - za[ip*dz];
                
              std::vector<double> cellprt_force = newton_force_(dens[i2]*cell_vol, disp); 

              axa[ip*dax] += cellprt_force[0];
              aya[ip*day] += cellprt_force[1];
              aza[ip*daz] += cellprt_force[2];

              tot_cell_force[0] += cellprt_force[0];
              tot_cell_force[1] += cellprt_force[1];
              tot_cell_force[2] += cellprt_force[2];

            }
          }
        }

        CkPrintf("Cell force (prt): %f, %f, %f\n", tot_cell_force[0], tot_cell_force[1], tot_cell_force[2]);

        std::vector<double> tot_prt_force (3, 0); // for debugging purposes

        // compute force sourced by other particles
        for (int ipt2 = 0; ipt2 < num_is_grav; ipt2++) {
          const int it2 = particle.type_index(particle_groups->item("is_gravitating",ipt2));

          int imass2 = 0;

          for (int ib2 = 0; ib2 < particle.num_batches(it2); ib2++) {

            const int np2 = particle.num_particles(it2,ib2);

            if (particle.has_attribute(it2,"mass")) {
              imass2 = particle.attribute_index(it2,"mass");
              prtmass2 = (enzo_float *) particle.attribute_array( it2, imass2, ib2);
              dm2 = particle.stride(it2,imass2);
            } 
            else {
              imass2 = particle.constant_index(it2,"mass");
              prtmass2 = (enzo_float*)particle.constant_value(it2,imass2);
              dm2 = 0;
            }

            const int ia_x2  = particle.attribute_index(it2,"x");
            const int ia_y2  = particle.attribute_index(it2,"y");
            const int ia_z2  = particle.attribute_index(it2,"z");

            enzo_float * xa2 =  (enzo_float *)particle.attribute_array (it2,ia_x2,ib2);
            enzo_float * ya2 =  (enzo_float *)particle.attribute_array (it2,ia_y2,ib2);
            enzo_float * za2 =  (enzo_float *)particle.attribute_array (it2,ia_z2,ib2);

            const int dx2 =  particle.stride(it2,ia_x2);
            const int dy2 =  particle.stride(it2,ia_y2);
            const int dz2 =  particle.stride(it2,ia_z2);

            for (int ip2=0; ip2 < np2; ip2++) {

              if (!(it == it2 && ib == ib2 && ip == ip2)) {
              
                // disp points from current particle to interacting particle
                std::vector<double> disp (3, 0);
                disp[0] = xa2[ip2*dx2] - xa[ip*dx]; 
                disp[1] = ya2[ip2*dy2] - ya[ip*dy];
                disp[2] = za2[ip2*dz2] - za[ip*dz];
                  
                std::vector<double> prtprt_force = newton_force_(prtmass2[ip2*dm2], disp); 

                axa[ip*dax] += prtprt_force[0];
                aya[ip*day] += prtprt_force[1];
                aza[ip*daz] += prtprt_force[2];

                tot_prt_force[0] += prtprt_force[0];
                tot_prt_force[1] += prtprt_force[1];
                tot_prt_force[2] += prtprt_force[2];
              }

            }
          }
        }

        CkPrintf("Particle force (prt): %f, %f, %f\n", tot_prt_force[0], tot_prt_force[1], tot_prt_force[2]);
        CkPrintf("Accel (prt): %f, %f, %f\n\n", axa[ip*dax], aya[ip*day], aza[ip*daz]);

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

void EnzoBlock::p_method_multipole_interact_direct (int fldsize, int prtsize, char * fldbuffer, char * prtbuffer)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->interact_direct_ (this, fldbuffer, prtbuffer);
}

void EnzoBlock::p_method_multipole_interact_direct_send (Index receiver)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  // method->interact_direct_send (this, receiver);
  method->pack_dens_ (this, receiver);
}

// // I don't think this is necessary? Can just call pack_dens_ in p_method_multipole_interact_direct_send
// void EnzoMethodMultipole::interact_direct_send(EnzoBlock * enzo_block, Index receiver) throw()
// {
//   pack_dens_(enzo_block, receiver);
// }


void EnzoMethodMultipole::pack_dens_(EnzoBlock * enzo_block, Index index_b) throw()
{
  Data * data = enzo_block->data();
  Field field = data->field();
  Particle particle = data->particle();

  enzo_float * dens = (enzo_float*) field.values ("density");

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;
  double lo[3];

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  enzo_block->cell_width(&hx, &hy, &hz);
  // hz = 1.0; // this is to make math easier -- delete later

  enzo_block->lower(lo, lo+1, lo+2);

  int fldsize = 0;

  SIZE_SCALAR_TYPE(fldsize, int, mx);
  SIZE_SCALAR_TYPE(fldsize, int, my);
  SIZE_SCALAR_TYPE(fldsize, int, mz);
  SIZE_SCALAR_TYPE(fldsize, double, hx);
  SIZE_SCALAR_TYPE(fldsize, double, hy);
  SIZE_SCALAR_TYPE(fldsize, double, hz);
  SIZE_ARRAY_TYPE(fldsize, double, lo, 3);
  SIZE_ARRAY_TYPE(fldsize, enzo_float, dens, mx*my*mz);

  char * fldbuffer = new char[fldsize];

  char * pc; 
  pc = fldbuffer;

  SAVE_SCALAR_TYPE(pc, int, mx);
  SAVE_SCALAR_TYPE(pc, int, my);
  SAVE_SCALAR_TYPE(pc, int, mz);
  SAVE_SCALAR_TYPE(pc, double, hx);
  SAVE_SCALAR_TYPE(pc, double, hy);
  SAVE_SCALAR_TYPE(pc, double, hz);
  SAVE_ARRAY_TYPE(pc, double, lo, 3);
  SAVE_ARRAY_TYPE(pc, enzo_float, dens, mx*my*mz);


  int prtsize = particle.data_size();
  char * prtbuffer = new char[prtsize];
  char * prtbuffer_next = particle.save_data(prtbuffer);

  enzo::block_array()[index_b].p_method_multipole_interact_direct (fldsize, prtsize, fldbuffer, prtbuffer);
}


void EnzoMethodMultipole::interact_direct_(Block * block, char * fldbuffer_b, char * prtbuffer_b) throw()
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
  // CkPrintf("ghosts: %d, %d, %d\n\n", gx, gy, gz);
  // CkPrintf("field dims: %d, %d, %d\n\n", mx, my, mz);

  block->cell_width(&hx, &hy, &hz);
  // hz = 1.0; // this is to make math easier -- delete later
  double cell_vol = hx * hy * hz;

  double lo[3];
  block->lower(lo, lo+1, lo+2);

  // unpacking data from fldbuffer_b
  int mx2, my2, mz2;
  double hx2, hy2, hz2;
  double lo2[3];

  char * pc;
  pc = fldbuffer_b;

  LOAD_SCALAR_TYPE(pc, int, mx2);
  LOAD_SCALAR_TYPE(pc, int, my2);
  LOAD_SCALAR_TYPE(pc, int, mz2);
  LOAD_SCALAR_TYPE(pc, double, hx2);
  LOAD_SCALAR_TYPE(pc, double, hy2);
  LOAD_SCALAR_TYPE(pc, double, hz2);
  LOAD_ARRAY_TYPE(pc, double, lo2, 3);

  double dens[mx2*my2*mz2];
  LOAD_ARRAY_TYPE(pc, double, dens, mx2*my2*mz2);


  Particle particle = data->particle();
  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();
  const int num_is_grav = particle_groups->size("is_gravitating");
  const int num_prtl_types = particle_descr->num_types();

  enzo_float * prtmass2 = NULL;
  int dm2;

  // unpacking particles from prtbuffer_b
  ParticleData new_p_data;
  Particle particle2 (particle_descr, &new_p_data);
  char * prtbuffer_next = particle2.load_data(prtbuffer_b);


  // loop over all cells in this Block
  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {

        int i = ix + mx*(iy + my*iz);

        CkPrintf("position: %f, %f, %f\n", (lo[0] + (ix-gx + 0.5)*hx), (lo[1] + (iy-gy + 0.5)*hy), (lo[2] + (iz-gz + 0.5)*hz));

        std::vector<double> tot_cell_force (3, 0); // for debugging

        // compute force sourced from cells in Block b
        for (int iz2 = gz; iz2 < mz2-gz; iz2++) {
          for (int iy2 = gy; iy2 < my2-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx2-gx; ix2++) {

              int i2 = ix2 + mx2*(iy2 + my2*iz2);

              // CkPrintf("position2: %f, %f, %f\n", (lo2[0] + (ix2-gx + 0.5)*hx2), (lo2[1] + (iy2-gy + 0.5)*hy2), (lo2[2] + (iz2-gz + 0.5)*hz2));

              // disp points from cell in current Block to cell in Block b
              std::vector<double> disp (3, 0);
              disp[0] = lo2[0] - lo[0] + (ix2-gx + 0.5)*hx2 - (ix-gx + 0.5)*hx; 
              disp[1] = lo2[1] - lo[1] + (iy2-gy + 0.5)*hy2 - (iy-gy + 0.5)*hy;
              disp[2] = lo2[2] - lo[2] + (iz2-gz + 0.5)*hz2 - (iz-gz + 0.5)*hz;

              // CkPrintf("disp: %f, %f, %f\n", disp[0], disp[1], disp[2]);
                
              std::vector<double> cell_force = newton_force_(dens[i2]*hx2*hy2*hz2, disp); 

              accel_x[i] += cell_force[0];
              accel_y[i] += cell_force[1];
              accel_z[i] += cell_force[2];

              tot_cell_force[0] += cell_force[0];
              tot_cell_force[1] += cell_force[1];
              tot_cell_force[2] += cell_force[2];

              // if (accel_z[i] != 0) {
              //   CkPrintf("position: %f, %f, %f\n", (lo[0] + (ix-gx + 0.5)*hx), (lo[1] + (iy-gy + 0.5)*hy), (lo[2] + (iz-gz + 0.5)*hz));
              //   CkPrintf("position2: %f, %f, %f\n", (lo2[0] + (ix2-gx + 0.5)*hx2), (lo2[1] + (iy2-gy + 0.5)*hy2), (lo2[2] + (iz2-gz + 0.5)*hz2));
              //   CkPrintf("cell force: %f, %f, %f\n", cell_force[0], cell_force[1], cell_force[2]);
              // }
            }
          }
        }
        CkPrintf("Leaf-leaf cell force: %f, %f, %f\n", tot_cell_force[0], tot_cell_force[1], tot_cell_force[2]);
        CkPrintf("Acceleration: %f, %f, %f\n\n", accel_x[i], accel_y[i], accel_z[i]);

        // compute cell force sourced from particles in Block b
        for (int ipt = 0; ipt < num_is_grav; ipt++) {
          // CkPrintf("are particles exerting force?\n");

          const int it = particle2.type_index(particle_groups->item("is_gravitating",ipt));

          int imass2 = 0;

          for (int ib = 0; ib < particle2.num_batches(it); ib++) {

            const int np = particle2.num_particles(it,ib);

            if (particle2.has_attribute(it,"mass")) {
              imass2 = particle2.attribute_index(it,"mass");
              prtmass2 = (enzo_float *) particle2.attribute_array( it, imass2, ib);
              dm2 = particle2.stride(it,imass2);
            } 
            else {
              imass2 = particle2.constant_index(it,"mass");
              prtmass2 = (enzo_float*)particle2.constant_value(it,imass2);
              dm2 = 0;
            }

            const int ia_x2  = particle2.attribute_index(it,"x");
            const int ia_y2  = particle2.attribute_index(it,"y");
            const int ia_z2  = particle2.attribute_index(it,"z");

            enzo_float * xa2 =  (enzo_float *)particle2.attribute_array (it,ia_x2,ib);
            enzo_float * ya2 =  (enzo_float *)particle2.attribute_array (it,ia_y2,ib);
            enzo_float * za2 =  (enzo_float *)particle2.attribute_array (it,ia_z2,ib);

            const int dx2 =  particle2.stride(it,ia_x2);
            const int dy2 =  particle2.stride(it,ia_y2);
            const int dz2 =  particle2.stride(it,ia_z2);

            for (int ip=0; ip < np; ip++) {

              // CkPrintf("This shouldn't print\n");
              
              // disp points from cell in current Block to particle in Block b
              std::vector<double> disp (3, 0);
              disp[0] = xa2[ip*dx2] - (lo[0] + (ix-gx + 0.5)*hx); 
              disp[1] = ya2[ip*dy2] - (lo[1] + (iy-gy + 0.5)*hy);
              disp[2] = za2[ip*dz2] - (lo[2] + (iz-gz + 0.5)*hz);
                
              std::vector<double> prtcell_force = newton_force_(prtmass2[ip*dm2], disp); 

              accel_x[i] += prtcell_force[0];
              accel_y[i] += prtcell_force[1];
              accel_z[i] += prtcell_force[2];

            }
          }
        }
      }
    }
  }


  // loop over all particles in this Block
  for (int it = 0; it < num_prtl_types; it++) {

    for (int ib = 0; ib < particle.num_batches(it); ib++) {

      const int np = particle.num_particles(it,ib);

      const int ia_x  = particle.attribute_index(it,"x");
      const int ia_y  = particle.attribute_index(it,"y");
      const int ia_z  = particle.attribute_index(it,"z");
      const int ia_ax  = particle.attribute_index(it,"ax");
      const int ia_ay  = particle.attribute_index(it,"ay");
      const int ia_az  = particle.attribute_index(it,"az");

      enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
      enzo_float * ya =  (enzo_float *)particle.attribute_array (it,ia_y,ib);
      enzo_float * za =  (enzo_float *)particle.attribute_array (it,ia_z,ib);
      enzo_float * axa =  (enzo_float *)particle.attribute_array (it,ia_ax,ib);
      enzo_float * aya =  (enzo_float *)particle.attribute_array (it,ia_ay,ib);
      enzo_float * aza =  (enzo_float *)particle.attribute_array (it,ia_az,ib);

      const int dx =  particle.stride(it,ia_x);
      const int dy =  particle.stride(it,ia_y);
      const int dz =  particle.stride(it,ia_z);
      const int dax =  particle.stride(it,ia_ax);
      const int day =  particle.stride(it,ia_ay);
      const int daz =  particle.stride(it,ia_az);

      for (int ip=0; ip < np; ip++) {
        
        // compute force sourced by cells in Block b
        for (int iz2 = gz; iz2 < mz2-gz; iz2++) {
          for (int iy2 = gy; iy2 < my2-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx2-gx; ix2++) {

              int i2 = ix2 + mx2*(iy2 + my2*iz2);

              // disp points from particle in current Block to cell in Block b
              std::vector<double> disp (3, 0);
              disp[0] = (lo2[0] + (ix2-gx + 0.5)*hx2) - xa[ip*dx]; 
              disp[1] = (lo2[1] + (iy2-gy + 0.5)*hy2) - ya[ip*dy];
              disp[2] = (lo2[2] + (iz2-gz + 0.5)*hz2) - za[ip*dz];
                
              std::vector<double> cellprt_force = newton_force_(dens[i2]*hx2*hy2*hz2, disp); 

              axa[ip*dax] += cellprt_force[0];
              aya[ip*day] += cellprt_force[1];
              aza[ip*daz] += cellprt_force[2];

            }
          }
        }

        // compute force sourced by particles in Block b
        for (int ipt2 = 0; ipt2 < num_is_grav; ipt2++) {
          const int it2 = particle2.type_index(particle_groups->item("is_gravitating",ipt2));

          int imass2 = 0;

          for (int ib2 = 0; ib2 < particle2.num_batches(it2); ib2++) {

            const int np2 = particle2.num_particles(it2,ib2);

            if (particle2.has_attribute(it2,"mass")) {
              imass2 = particle2.attribute_index(it2,"mass");
              prtmass2 = (enzo_float *) particle2.attribute_array( it2, imass2, ib2);
              dm2 = particle2.stride(it2,imass2);
            } 
            else {
              imass2 = particle2.constant_index(it2,"mass");
              prtmass2 = (enzo_float*)particle2.constant_value(it2,imass2);
              dm2 = 0;
            }

            const int ia_x2  = particle2.attribute_index(it2,"x");
            const int ia_y2  = particle2.attribute_index(it2,"y");
            const int ia_z2  = particle2.attribute_index(it2,"z");

            enzo_float * xa2 =  (enzo_float *)particle2.attribute_array (it2,ia_x2,ib2);
            enzo_float * ya2 =  (enzo_float *)particle2.attribute_array (it2,ia_y2,ib2);
            enzo_float * za2 =  (enzo_float *)particle2.attribute_array (it2,ia_z2,ib2);

            const int dx2 =  particle2.stride(it2,ia_x2);
            const int dy2 =  particle2.stride(it2,ia_y2);
            const int dz2 =  particle2.stride(it2,ia_z2);

            for (int ip2=0; ip2 < np2; ip2++) {
              
              // disp points from particle in current Block to particle in Block b
              std::vector<double> disp (3, 0);
              disp[0] = xa2[ip2*dx2] - xa[ip*dx]; 
              disp[1] = ya2[ip2*dy2] - ya[ip*dy];
              disp[2] = za2[ip2*dz2] - za[ip*dz];
                
              std::vector<double> prtprt_force = newton_force_(prtmass2[ip2*dm2], disp); 

              axa[ip*dax] += prtprt_force[0];
              aya[ip*day] += prtprt_force[1];
              aza[ip*daz] += prtprt_force[2];

            }
          }
        }

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

/**** for particle debugging  ****/

// see EnzoInitialIsolatedGalaxy for more detailed initialization

 void EnzoMethodMultipole::InitializeParticles(Block * block, int nprtls, double prtls[][7]){

  Particle particle = block->data()->particle();
  ParticleDescr * particle_descr = cello::particle_descr();

  double lo[3];
  double hi[3];
  block->lower(lo, lo+1, lo+2);
  block->upper(hi, hi+1, hi+2);

  // CkPrintf("lo: %f, %f, %f\n", lo[0], lo[1], lo[2]);
  // CkPrintf("hi: %f, %f, %f\n", hi[0], hi[1], hi[2]);

  //
  // Loop through all particle types and initialize their positions and
  // accelerations.
  //

  int ntypes = 2;
  int * particleIcTypes = new int[ntypes];

  particleIcTypes[0] = particle_descr->type_index("star");
  particleIcTypes[1] = particle_descr->type_index("trace");

  // particleIcFileNames.push_back("halo.dat");
  // nparticles_ = std::max(nparticles_, nlines("halo.dat"));
  // ipt++;

  // Loop over all particle types and initialize
  for(int ipt = 0; ipt < ntypes; ipt++){

    int it   = particleIcTypes[ipt];

    // obtain particle attribute indexes for this type
    int ia_x = particle.attribute_index (it, "x");
    int ia_y = particle.attribute_index (it, "y");
    int ia_z = particle.attribute_index (it, "z");
    int ia_ax = particle.attribute_index (it, "ax");
    int ia_ay = particle.attribute_index (it, "ay");
    int ia_az = particle.attribute_index (it, "az");

    int ib  = 0; // batch counter
    int ipp = 0; // particle counter

    // this will point to the particular value in the
    // particle attribute array

    enzo_float * px   = 0;
    enzo_float * py   = 0;
    enzo_float * pz   = 0;
    enzo_float * pax  = 0;
    enzo_float * pay  = 0;
    enzo_float * paz  = 0;

    if (ipt != ntypes-1) {

      int ia_m = particle.attribute_index (it, "mass");
      enzo_float * prtmass = 0;

      for (int i = 0; i < nprtls-1; i++){

        // CkPrintf("pos: %f, %f, %f\n", prtls[i][1], prtls[i][2], prtls[i][3]);

        if (prtls[i][1] >= lo[0] && prtls[i][1] < hi[0] &&
            prtls[i][2] >= lo[1] && prtls[i][2] < hi[1] &&
            prtls[i][3] >= lo[2] && prtls[i][3] < hi[2]) {

          // CkPrintf("added a particle!\n");

          int new_particle = particle.insert_particles(it, 1);
          particle.index(new_particle,&ib,&ipp);

          // get pointers to each of the associated arrays
          prtmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
          px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
          py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
          pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
          pax   = (enzo_float *) particle.attribute_array(it, ia_ax, ib);
          pay   = (enzo_float *) particle.attribute_array(it, ia_ay, ib);
          paz   = (enzo_float *) particle.attribute_array(it, ia_az, ib);


          // set the particle values
          prtmass[ipp] = prtls[i][0];
          px[ipp]    = prtls[i][1];
          py[ipp]    = prtls[i][2];
          pz[ipp]    = prtls[i][3];
          pax[ipp]   = prtls[i][4];
          pay[ipp]   = prtls[i][5];
          paz[ipp]   = prtls[i][6];
        }
      }
    } 

    else {

      int i = nprtls - 1;

      if (prtls[i][1] >= lo[0] && prtls[i][1] < hi[0] &&
            prtls[i][2] >= lo[1] && prtls[i][2] < hi[1] &&
            prtls[i][3] >= lo[2] && prtls[i][3] < hi[2]) {

          // CkPrintf("added a particle!\n");

          int new_particle = particle.insert_particles(it, 1);
          particle.index(new_particle,&ib,&ipp);

          // get pointers to each of the associated arrays
          px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
          py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
          pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
          pax   = (enzo_float *) particle.attribute_array(it, ia_ax, ib);
          pay   = (enzo_float *) particle.attribute_array(it, ia_ay, ib);
          paz   = (enzo_float *) particle.attribute_array(it, ia_az, ib);


          // set the particle values
          px[ipp]    = prtls[i][1];
          py[ipp]    = prtls[i][2];
          pz[ipp]    = prtls[i][3];
          pax[ipp]   = prtls[i][4];
          pay[ipp]   = prtls[i][5];
          paz[ipp]   = prtls[i][6];
      }
    }

  } // end loop over particle types


  return;
}