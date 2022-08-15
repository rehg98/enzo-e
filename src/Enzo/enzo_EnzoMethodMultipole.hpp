// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMultipole.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu) 
/// @date     Tuesday July 19, 2022
/// @brief    [\ref Enzo] Declaration of EnzoMethodMultipole
///           Compute multipoles and pass multipoles up octree

#ifndef ENZO_ENZO_METHOD_MULTIPOLE_HPP
#define ENZO_ENZO_METHOD_MULTIPOLE_HPP

class EnzoMethodMultipole : public Method {

  /// @class    EnzoMethodMultipole
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Compute multipoles and pass multipoles up octree

public: // interface

  /// Create a new EnzoMethodMultipole object
  EnzoMethodMultipole(double timeStep, int maxLevel);

  EnzoMethodMultipole()
    : Method(),
      timeStep_(10000.0),
      maxLevel_(5)
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodMultipole);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodMultipole (CkMigrateMessage *m)
    : Method (m),
      timeStep_(10000.0),
      maxLevel_(5),
      i_sync_restrict_(-1),
      i_msg_restrict_()
  {for (int i = 0; i < cello::num_children(); i++) i_msg_restrict_[i] = -1; }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
     
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "multipole"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();

  void restrict_send (EnzoBlock * enzo_block) throw();

  void restrict_recv (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();

  MultipoleMsg * pack_multipole_ (EnzoBlock * enzo_block) throw();

  void unpack_multipole_ (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();

 double * pmass(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr,i_mass_);
  }


  /// returns an array of three elements corresponding to components of COM
  double * pcom(Block * block)
  {
    
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr,i_com_);
  }

  /// returns a 9-element array corresponding to the components of the quadrupole tensor
  double * pquadrupole(Block * block)
  {
    
    
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr,i_quadrupole_);
  }



  /// Access the restrict Sync Scalar value for the Block
  Sync * psync_restrict(Block * block)
  {
    ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
    ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
    return scalar_data->value(scalar_descr,i_sync_restrict_);
  }

  /// Access the Field message for buffering restriction data
  MultipoleMsg ** pmsg_restrict(Block * block, int ic)
  {
    
    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();
    ScalarDescr *        scalar_descr = cello::scalar_descr_void();

    return (MultipoleMsg **)scalar_data->value(scalar_descr,i_msg_restrict_[ic]);
  }

protected: // methods

  void compute_ (Block * block) throw();

  void compute_multipoles_ (Block * block) throw();

  void begin_cycle_ (EnzoBlock * enzo_block) throw();

  double * shift_quadrupole_(double *oldQuadrupole, double totMass, double *oldCOM, double *newCOM) throw();


protected: // attributes

  /// Maximum timestep
  double timeStep_;

  /// Maximum level
  int maxLevel_;

  /// multipoles for this block
  int i_mass_;
  int i_com_;
  int i_quadrupole_;

  int i_sync_restrict_;
  int i_msg_restrict_[8];
};

#endif /* ENZO_ENZO_METHOD_MULTIPOLE_HPP */
