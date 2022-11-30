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

public: // interface -- which methods should be public and which protected?

  /// Create a new EnzoMethodMultipole object
  EnzoMethodMultipole(double timeStep, double theta);

  EnzoMethodMultipole()
    : Method(),
      timeStep_(10000.0),
      theta_(0),
      is_volume_(-1),
      block_volume_(),
      max_volume_(0)
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodMultipole);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodMultipole (CkMigrateMessage *m)
    : Method (m),
      timeStep_(10000.0),
      i_sync_restrict_(-1),
      //i_sync_prolong_(-1),
      i_msg_restrict_(),
      i_msg_prolong_(-1),
      theta_(0),
      is_volume_(-1),
      block_volume_(),
      max_volume_(0)
  { for (int i = 0; i < cello::num_children(); i++) i_msg_restrict_[i] = -1; }



  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
     
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "multipole"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();

  void begin_up_cycle_ (EnzoBlock * enzo_block) throw();
  void begin_down_cycle_ (EnzoBlock * enzo_block) throw();

  void restrict_send (EnzoBlock * enzo_block) throw();
  void prolong_send (EnzoBlock * enzo_block) throw();

  void restrict_recv (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();
  void prolong_recv (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();

  MultipoleMsg * pack_multipole_ (EnzoBlock * enzo_block) throw();
  MultipoleMsg * pack_coeffs_ (EnzoBlock * enzo_block) throw();

  void unpack_multipole_ (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();
  void unpack_coeffs_ (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();

  /*****************************************************/

  /// Access the restrict Sync Scalar value for the Block
  Sync * psync_restrict(Block * block)
  {
    ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
    ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
    return scalar_data->value(scalar_descr, i_sync_restrict_);
  }

  /// Access the prolong Sync Scalar value for the Block
  //  Sync * psync_prolong(Block * block)
  // {
  //  ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
  //  ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
  //  return scalar_data->value(scalar_descr, i_sync_prolong_);
  //}

  /// Access the Field message for buffering restriction data
  MultipoleMsg ** pmsg_restrict(Block * block, int ic)
  {
    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();
    ScalarDescr *        scalar_descr = cello::scalar_descr_void();
    return (MultipoleMsg **)scalar_data->value(scalar_descr, i_msg_restrict_[ic]);
  }

  /// Access the Field message for buffering prolongation data
  MultipoleMsg ** pmsg_prolong(Block * block)
  {
    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();
    ScalarDescr *        scalar_descr = cello::scalar_descr_void();
    
    void ** data = scalar_data->value(scalar_descr, i_msg_prolong_);

    if (data == nullptr){
      ERROR("data is null", "data is null");
    }

    return (MultipoleMsg **)data;
  }

  double * pmass(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_mass_);
  }

  /// returns an array of three elements corresponding to components of COM
  double * pcom(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_com_);
  }

  /// returns a 9-element array corresponding to the components of the quadrupole tensor
  double * pquadrupole(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_quadrupole_);
  }

  /// returns a 3-element array
  double * pc1(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_c1_);
  }

  /// returns a 9-element array
  double * pc2(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_c2_);
  }

  /// returns a 27-element array
  double * pc3(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_c3_);
  }

  /*****************************************************/


protected: // methods

  void compute_ (Block * block) throw();


  /**********    FMM specific functions       **********/

  void compute_multipoles_ (Block * block) throw();

  void evaluate_force_ (Block * block) throw();  // compute force from taylor coeffs

  // compute the Newtonian acceleration induced on this cell/particle by Cell/Particle B
  std::vector<double> newton_force_ (double mass_b, std::vector<double> displacement) throw()
  {
    std::vector<double> accel_vec (3, 0);

    double disp_norm = sqrt(dot_11_(displacement, displacement));

    double accel_mag = mass_b / (disp_norm * disp_norm);   // how does the code represent G?
    accel_vec = dot_scalar_(accel_mag / disp_norm, displacement, 3);

    return accel_vec;
  }


  /********** functions for dual tree walk **********/

public:

  void interact_approx_ (Block * block, MultipoleMsg * msg_b) throw(); // compute Taylor coeffs for two interacting blocks
  void interact_direct_ (Block * block, MultipoleMsg * msg_b) throw(); // compute Taylor coeffs for two interacting blocks

  void interact_approx_send(EnzoBlock * enzo_block, Index receiver) throw();
  void interact_direct_send(EnzoBlock * enzo_block, Index receiver) throw();

  void dual_walk_ (EnzoBlock * enzo_block) throw();

  void traverse (EnzoBlock * enzo_block, Index index_b, int type_b);

  /// Identified pair of blocks considered to be "distant"
  void traverse_approx_pair (EnzoBlock * enzo_block,
                            Index index_a, int volume_a,
                            Index index_b, int volume_b);
                              
  /// Identified pair of blocks considered to be "close". (Block
  /// A may be the same as block B but won't be called twice.)
  void traverse_direct_pair (EnzoBlock * enzo_block,
                            Index index_a, int volume_a,
                            Index index_b, int volume_b);

  MultipoleMsg * pack_dens_ (EnzoBlock * enzo_block) throw(); 

  void update_volume (EnzoBlock * enzo_block, Index index, int volume);

protected:

  bool is_far_ (EnzoBlock * enzo_block, Index index_b,
              double * ra, double * rb) const;

  long long * volume_(Block * block) {
    return block->data()->scalar_long_long().value(is_volume_);
  }


  /**********   convenience functions for tensor arithmetic ***********/

  std::vector<double> shift_quadrupole_(double * old_quadrupole, double tot_mass, double * old_com, double * new_com) throw()
  {
    std::vector<double> new_com_vec (new_com, new_com + 3);
    std::vector<double> old_com_vec (old_com, old_com + 3);

    std::vector<double> disp = subtract_(new_com_vec, old_com_vec, 3);
    
    std::vector<double> new_quadrupole (9, 0);
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        new_quadrupole[3*i + j] = old_quadrupole[3*i + j] + tot_mass * disp[i] * disp[j];
      }
    }
    
    return new_quadrupole;

  }

  // subtract two arrays
  std::vector<double> subtract_(std::vector<double> a, std::vector<double> b, int size) throw()
  {
    // should assert that a and b have size "size"
    std::vector<double> diff (size, 0);
    
    for (int i = 0; i < size; i++) {
      diff[i] = a[i] - b[i];
    }

    return diff;
    
  }

  // add two arrays
  std::vector<double> add_(std::vector<double> a, std::vector<double> b, int size) throw()
  {
    // should assert that a and b have size "size"
    std::vector<double> sum (size, 0);

    for (int i = 0; i < size; i++) {
      sum[i] = a[i] + b[i];
    }

    return sum;
  }

  // compute the outer product of two (3-elements) vectors 
  std::vector<double> outer_(std::vector<double> a, std::vector<double> b) throw()
  {
    // assert that a and b are 3-element vectors 
    // double * prod = new double[9];

    std::vector<double> prod (9, 0);

    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        prod[3*i + j] = a[i] * b[j];
      }
    }

    return prod;
  }

  // multiply an array by a scalar
  std::vector<double> dot_scalar_(double a, std::vector<double> b, int size) throw()
  {
    // double * prod = new double[size];
    std::vector<double> prod (size, 0);

    for (int i = 0; i < size; i++) {
      prod[i] = a * b[i];
    }

    return prod;
  }

  // compute the dot product of two vectors
  double dot_11_(std::vector<double> a, std::vector<double> b) throw()
  {
    double prod = 0;

    for (int i = 0; i < 3; i++) {
      prod += a[i] * b[i];
    }

    return prod;
  }

  // compute the dot product of a vector 'a' with a rank-2 tensor 'b'
  std::vector<double> dot_12_(std::vector<double> a, std::vector<double> b) throw()
  {
    // double * prod = new double[3];
    std::vector<double> prod (3, 0);

    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        prod[j] += a[i] * b[3*i + j];
      }
    }

    return prod;
  }

  // compute the dot product of a vector 'a' with a rank-3 tensor 'b'
  std::vector<double> dot_13_(std::vector<double> a, std::vector<double> b) throw()
  {

    // double * prod = new double[9];
    std::vector<double> prod (9, 0);

    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          prod[3*j + k] += a[i] * b[9*k + 3*i + j];
        }
      }
    }

    return prod;

  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-3 tensor 'b'
  std::vector<double> dot_23_(std::vector<double> a, std::vector<double> b) throw()
  {
    // double * prod = new double[3];
    std::vector<double> prod (3, 0);

    for (int k = 0; k < 3; k++) {
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          prod[k] += a[3*i + j] * b[9*k + 3*i + j];
        }
      }
    }

    return prod;
  }


protected: // attributes

  /// Maximum timestep
  double timeStep_;
  
  // int i_sync_prolong_;
  int i_sync_restrict_;
  int i_msg_prolong_;
  int i_msg_restrict_[8];

  /// indices to multipoles for this block
  int i_mass_;
  int i_com_;
  int i_quadrupole_;

  /// indices to Taylor coefficients for this block
  int i_c1_;
  int i_c2_;
  int i_c3_;

  /// Parameter controlling the multipole acceptance criterion, defined
  /// as theta * distance(A,B) > (radius(A) + radius(B))
  double theta_;

  // /// Minimum/maximum mesh refinement level (saved for efficiency)
  // int min_level_;
  // int max_level_;

  /// ScalarData<long long> index for summing block volumes to
  /// test for termination
  int is_volume_;

  /// volume of block in each level, with finest level normalized to 1. Used
  /// to test for termination
  std::vector<int> block_volume_;

  /// Volume of domain in terms of finest blocks
  long long max_volume_;

};

#endif /* ENZO_ENZO_METHOD_MULTIPOLE_HPP */
