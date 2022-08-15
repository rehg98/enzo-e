// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMultipole.cpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     July 19, 2022
/// @brief    Compute multipoles and pass multipoles up octree

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodMultipole::EnzoMethodMultipole (double timeStep, int maxLevel)
  : Method(),
    timeStep_(timeStep),
    maxLevel_(maxLevel)  
{ 
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  i_sync_restrict_ = scalar_descr_sync->new_value("multipole:restrict");

  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  for (int ic=0; ic<cello::num_children(); ic++) {
    i_msg_restrict_[ic] = scalar_descr_void->new_value("multipole:msg_restrict");
  }

  cello::define_field ("density");

  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("density");

  ScalarDescr * scalar_descr_double = cello::scalar_descr_double();
  i_mass_ = scalar_descr_double->new_value("multipole:mass");

  i_com_ = scalar_descr_double->new_value("multipole:com", 3);
    
  i_quadrupole_ = scalar_descr_double->new_value("multipole:quadrupole", 9);


    

  
}

//----------------------------------------------------------------------

void EnzoMethodMultipole::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | timeStep_;
  p | maxLevel_;
  p | i_sync_restrict_;
  p | i_mass_;
  p | i_com_;
  p | i_quadrupole_;
  PUParray(p, i_msg_restrict_, 8);
  
}

//----------------------------------------------------------------------

void EnzoMethodMultipole::compute ( Block * block) throw()
{

  Sync * sync_restrict = psync_restrict(block);
  sync_restrict->set_stop(1 + cello::num_children());
  
  double * mass = pmass(block);
  double * com = pcom(block);
  double * quadrupole = pquadrupole(block);

  
  *mass = 0;

  for (int i = 0; i < 3; i++){
    com[i]  = 0;
  }
  
  for (int i = 0; i < 9; i++){
    quadrupole[i]  = 0;
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

void EnzoMethodMultipole::compute_multipoles_ (Block * block) throw()
{
  Data * data = block->data();
  Field field = data->field();

  EnzoBlock * enzo_block = enzo::block(block);
  int level = enzo_block->level();
  // CkPrintf("multipole level: %d\n", level);

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

  double * mass_ = pmass(block);
  double * com_ = pcom(block);
  double * quadrupole_ = pquadrupole(block);

  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {
        
        double dens = density[ix + mx * (iy + iz * my)];
        // CkPrintf("dens: %f\n", dens);

        *mass_ += dens;

        // CkPrintf("x coord: %f\n", (lo[0] + (ix-gx + 0.5)*hx));
        // CkPrintf("y coord: %f\n", (lo[1] + (iy-gy + 0.5)*hy));

	      weighted_sum[0] += dens * (lo[0] + (ix-gx + 0.5)*hx); 
	      weighted_sum[1] += dens * (lo[1] + (iy-gy + 0.5)*hy); 
	      weighted_sum[2] += dens * (lo[2] + (iz-gz + 0.5)*hz);

      }
    }
  }

  com_[0] = weighted_sum[0] / *mass_;
  com_[1] = weighted_sum[1] / *mass_;
  com_[2] = weighted_sum[2] / *mass_;

  *mass_ *= cell_vol;
  // CkPrintf("Block mass: %f\n", mass_);
  

  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {
        
        double dens = density[ix + mx * (iy + iz * my)];
        //CkPrintf("dens: %f\n", dens);

        double disp[3];
        disp[0] = (lo[0] + (ix-gx + 0.5)*hx) - com_[0]; 
        disp[1] = (lo[1] + (iy-gy + 0.5)*hy) - com_[1];
        disp[2] = (lo[2] + (iz-gz + 0.5)*hz) - com_[2]; 

        //CkPrintf("disp[1]: %f\n", disp[1]);

        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            quadrupole_[3*i + j] += dens * cell_vol * disp[i] * disp[j];
          }
        }

      }
    }
  }


}



void EnzoMethodMultipole::compute_ (Block * block) throw()
{
  EnzoBlock * enzo_block = enzo::block(block);
  Sync * sync_restrict = psync_restrict(block);

  double * mass_ = pmass(block);
  double * com_ = pcom(block);
  double * quadrupole_ = pquadrupole(block);

  int level = enzo_block->level();
  
  if (block->is_leaf()) {
    
    compute_multipoles_ (block);
    begin_cycle_ (enzo_block);
    
  } else {

    if (0 <= level && level < maxLevel_) {
      CkPrintf("Non-leaf mass: %f\n", *mass_);
      restrict_recv(enzo_block, nullptr);
    }
  }
    
}


void EnzoMethodMultipole::begin_cycle_(EnzoBlock * enzo_block) throw()
{
  const int level = enzo_block->level();

  double * mass_ = pmass(enzo_block);
  double * com_ = pcom(enzo_block);
  double * quadrupole_ = pquadrupole(enzo_block);

  if (level == 0) {
    CkPrintf("total mass: %f\n", *mass_); 
    CkPrintf("COM: (%f, %f, %f)\n", com_[0], com_[1], com_[2]);
    CkPrintf("quadrupole: ");
    for (int i = 0; i < 9; i++) {
      CkPrintf("%f  ", quadrupole_[i]);
    }
    CkPrintf("\n");
  }
  else {
    // CkPrintf("cycle level: %d\n", level);
    restrict_send (enzo_block);
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
  double * mass_ = pmass(enzo_block);
  double * com_ = pcom(enzo_block);
  double * quadrupole_ = pquadrupole(enzo_block);
  
  int level = enzo_block->level();
  
  CkPrintf("recv start mass: %f\n", *mass_);
  CkPrintf("recv start level: %d\n", level);

  // Save field message from child
  if (msg != nullptr)
    {
      CkPrintf("recv pre-pmsg mass: %f\n", *mass_);
      *pmsg_restrict(enzo_block,msg->child_index()) = msg;
      CkPrintf("recv post-pmsg mass: %f\n", *mass_);
    }
  
  // Continue if all expected messages received
  if (psync_restrict(enzo_block)->next() ) {

    //CkPrintf("recv pre-loop mass: %f\n", *mass_);
    // Restore saved messages then clear
    for (int i=0; i<cello::num_children(); i++) {
      msg = *pmsg_restrict(enzo_block,i);
      *pmsg_restrict(enzo_block,i) = nullptr;

      // Unpack multipoles from message then delete message
      //CkPrintf("unpacking\n");
      unpack_multipole_(enzo_block,msg);
    }
  
    begin_cycle_ (enzo_block);
  }
  CkPrintf("recv end mass: %f\n", *mass_);

  
}


MultipoleMsg * EnzoMethodMultipole::pack_multipole_(EnzoBlock * enzo_block) throw()
{
 
  // Create a MultipoleMsg for sending data to parent
 
  MultipoleMsg * msg  = new MultipoleMsg;

  double * mass_ = pmass(enzo_block);
  double * com_ = pcom(enzo_block);
  double * quadrupole_ = pquadrupole(enzo_block);
 
  msg->mass = *mass_;

  for (int i = 0; i < 3; i++) {
    msg->com[i] = com_[i];
  }

  for (int i = 0; i < 9; i++) {
    msg->quadrupole[i] = quadrupole_[i];
  }

  const int level = enzo_block->level();
  if (level > 0) {
    enzo_block->index().child(level,msg->ic3,msg->ic3+1,msg->ic3+2);
  }

  return msg;

}

double * EnzoMethodMultipole::shift_quadrupole_(double *oldQuadrupole, double totMass, double *oldCOM, double *newCOM) throw()
{
  double disp[3];
  for (int i = 0; i < 3; i++) {
    disp[i] = newCOM[i] - oldCOM[i];
  }
   
  double * newQuadrupole = new double[9];
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      newQuadrupole[3*i + j] = oldQuadrupole[3*i + j] + totMass * disp[i] * disp[j];
    }
  }
   
  return newQuadrupole;
}

void EnzoMethodMultipole::unpack_multipole_
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{
  // copy data from msg to this EnzoBlock

  double * mass_ = pmass(enzo_block);
  double * com_ = pcom(enzo_block);
  double * quadrupole_ = pquadrupole(enzo_block);

  CkPrintf("Unpack start mass: %f\n", *mass_);
  
  double mass = msg->mass;
  double * com = msg->com;
  double * quadrupole = msg->quadrupole;

  double newCOM[3];
  for (int i = 0; i < 3; i++) {
    newCOM[i] = (*mass_ * com_[i] + mass * com[i]) / (*mass_ + mass);
  }
   
  double *shiftedQThis = shift_quadrupole_(quadrupole_, *mass_, com_, newCOM);
  double *shiftedQThat = shift_quadrupole_(quadrupole, mass, com, newCOM);
  for (int i = 0; i < 9; i++) {
      quadrupole_[i] = shiftedQThis[i] + shiftedQThat[i];
  }
  delete[] shiftedQThis;
  delete[] shiftedQThat;

  for (int i = 0; i < 3; i++)
    com_[i] = newCOM[i];

  
  *mass_ += mass;

  CkPrintf("Unpack end mass: %f\n", *mass_);
  delete msg;
}
