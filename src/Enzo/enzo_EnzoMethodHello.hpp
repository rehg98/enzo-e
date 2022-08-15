// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHello.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu) 
/// @date     Tuesday May 10, 2022
/// @brief    [\ref Enzo] Declaration of EnzoMethodHello
///           Simple Hello World Method

#ifndef ENZO_ENZO_METHOD_HELLO_HPP
#define ENZO_ENZO_METHOD_HELLO_HPP

class EnzoMethodHello : public Method {

  /// @class    EnzoMethodHello
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Simple Hello World method 

public: // interface

  /// Create a new EnzoMethodHello object
  EnzoMethodHello(double timeStep, int maxLevel);

  EnzoMethodHello()
    : Method(),
      timeStep_(10000.0),
      maxLevel_(5)
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodHello);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodHello (CkMigrateMessage *m)
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
  { return "hello"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();

  void restrict (EnzoBlock * enzo_block) throw();

  void restrict_send (EnzoBlock * enzo_block) throw();

  void restrict_recv (EnzoBlock * enzo_block, FieldMsg * msg) throw();

  FieldMsg * pack_residual_ (EnzoBlock * enzo_block) throw();

  void unpack_residual_ (EnzoBlock * enzo_block, FieldMsg * msg) throw();

  //void p_method_hello_restrict();

  //void p_method_hello_restrict_recv(FieldMsg * msg);

   /// Access the restrict Sync Scalar value for the Block
  Sync * psync_restrict(Block * block)
  {
    ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
    ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
    return scalar_data->value(scalar_descr,i_sync_restrict_);
  }

  /// Access the Field message for buffering restriction data
  FieldMsg ** pmsg_restrict(Block * block, int ic)
  {
    // check if scalar data is being set to null pointer
    //CkPrintf("pmsg_restrict: child index = %d\n", ic);
    //CkPrintf("pmsg i_msg = %d\n", i_msg_restrict_[0]);
    
    //CkPrintf("pmsg_restrict: scalar data = %d\n", block->data()->scalar_data_void());

    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();

    //CkPrintf("pmsg_restrict: scalar_data = %d\n", scalar_data);

    ScalarDescr *        scalar_descr = cello::scalar_descr_void();

    //CkPrintf("pmsg_restrict: scalar_descr = %d\n", scalar_descr);
    //CkPrintf("pmsg_restrict: i_msg_restrict = %d\n", i_msg_restrict_[ic]);

    return (FieldMsg **)scalar_data->value(scalar_descr,i_msg_restrict_[ic]);
  }

protected: // methods

  void compute_ (Block * block) throw();

  void begin_cycle_ (EnzoBlock * enzo_block) throw();

protected: // attributes

  /// Maximum timestep
  double timeStep_;

  /// Maximum level
  int maxLevel_;

  /// Number to pass between levels
  int n_;

  int i_sync_restrict_;
  int i_msg_restrict_[8];
};

#endif /* ENZO_ENZO_METHOD_HELLO_HPP */
