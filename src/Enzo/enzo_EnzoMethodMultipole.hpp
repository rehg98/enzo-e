// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMultipole.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu) 
/// @date     Tuesday July 19, 2022
/// @brief    [\ref Enzo] Declaration of EnzoMethodMultipole
/// Compute accelerations on gas and particles using FMM

#ifndef ENZO_ENZO_METHOD_MULTIPOLE_HPP
#define ENZO_ENZO_METHOD_MULTIPOLE_HPP

#include "EnzoMethodEwald.hpp"

class EnzoMethodMultipole : public Method {

  /// @class    EnzoMethodMultipole
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Compute multipoles and pass multipoles up octree

public: // interface -- which methods should be public and which protected?

  /// Create a new EnzoMethodMultipole object
  EnzoMethodMultipole(double timeStep, double theta, double eps0, double r0,
                      int interp_xpoints, int interp_ypoints, int interp_zpoints);

  EnzoMethodMultipole()
    : Method(),
      timeStep_(10000.0),
      theta_(0),
      eps0_(0),
      r0_(0),
      is_volume_(-1),
      block_volume_(),
      max_volume_(0),
      interp_xpoints_(16),
      interp_ypoints_(16),
      interp_zpoints_(16),
      ewald_(nullptr)
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
      eps0_(0),
      r0_(0),
      is_volume_(-1),
      block_volume_(),
      max_volume_(0),
      interp_xpoints_(16),
      interp_ypoints_(16),
      interp_zpoints_(16),
      ewald_(nullptr)
  { for (int i = 0; i < cello::num_children(); i++) i_msg_restrict_[i] = -1; }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  ~EnzoMethodMultipole()
  {
    if (ewald_ != nullptr) {
      delete ewald_;
    }
  }

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

  /// returns a 6-element array corresponding to the components of the quadrupole tensor
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

  /// returns a 6-element array
  double * pc2(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_c2_);
  }

  /// returns a 10-element array
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

  // // compute the Newtonian acceleration induced on this cell/particle by Cell/Particle B
  // std::vector<double> newton_force_ (double mass_b, std::vector<double> displacement) throw()
  // {
  //   std::vector<double> accel_vec (3, 0);

  //   double disp_norm = sqrt(dot_11_(displacement, displacement));

  //   double accel_mag = mass_b / (disp_norm * disp_norm);   // how does the code represent G?
  //   accel_vec = dot_scalar_(accel_mag / disp_norm, displacement, 3);

  //   return accel_vec;
  // }

  // compute the Newtonian acceleration induced on this cell/particle by Cell/Particle B;
  // force is softened with gravitational softening length eps0 and cutoff distance r0
  std::vector<double> newton_force_ (double mass_b, std::vector<double> displacement) throw()
  {
    std::vector<double> accel_vec (3, 0);

    double disp_norm = sqrt(dot_11_(displacement, displacement));  // norm of displacement vector
    double eps = epsilon(disp_norm, eps0_, r0_);                   // softening
    double soft_disp = disp_norm + eps;                            // softened displacement

    double accel_mag = mass_b / (soft_disp * soft_disp);   
    accel_vec = dot_scalar_(accel_mag / soft_disp, displacement, 3);

    return accel_vec;
  }

  /* Gravitational softening law from Springel et al. 2013; r is the displacement and eps0 is the softening length. 
      The softening has a finite range, going to 0 for distances greater than r0 (with r0 being smaller than half the smallest
      box dimension)  */
   double epsilon(double r, double eps0, double r0)
   {
    if (r >= r0)
      return 0; 
      
    else 
      return -2.8 * eps0 / W2(r / (2.8 * eps0)) - r; 
   }

  /* Kernel for gravitational softening (from Springel, Yoshida, White 2001) */
   double W2(double u)
   {
    if ((u >= 0) && (u < 0.5))
      return 16./3 * pow(u,2) - 48./5 * pow(u,4) + 32./5 * pow(u,5) - 14./5;
      
    else if ((u >= 0.5) && (u < 1))
      return 1./15 * pow(u,-1) + 32./3 * pow(u,2) - 16 * pow(u,3) + 48./5 * pow(u,4) - 32./15 * pow(u,5) - 16./5;
      
    else 
      return -1./u;
   }


  /********** functions for dual tree walk **********/

public:

  void interact_approx_ (Block * block, MultipoleMsg * msg_b) throw(); // compute Taylor coeffs for two interacting blocks
  void interact_direct_ (Block * block, char * fldbuffer_b, char * prtbuffer_b) throw(); // compute Newtonian force directly for two interacting leaves

  void interact_approx_send(EnzoBlock * enzo_block, Index receiver) throw();
  // void interact_direct_send(EnzoBlock * enzo_block, Index receiver) throw();

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

  void pack_dens_ (EnzoBlock * enzo_block, Index index_b) throw(); 

  void update_volume (EnzoBlock * enzo_block, Index index, int volume);

protected:

  bool is_far_ (EnzoBlock * enzo_block, Index index_b,
              double * ra, double * rb) const;

  long long * volume_(Block * block) {
    return block->data()->scalar_long_long().value(is_volume_);
  }

  /**********    functions for debugging particles          ************/
  void InitializeParticles (Block * block, int nprtls, double prtls[][7]);


public:
  /**********   convenience functions for tensor arithmetic ***********/

  // to what extent should i consider replacing std::vectors with std::array? Don't know size for some methods

  std::array<double, 6> shift_quadrupole_(double * old_quadrupole, double tot_mass, double * old_com, double * new_com) throw()
  {
    // std::vector<double> new_com_vec (new_com, new_com + 3);
    // std::vector<double> old_com_vec (old_com, old_com + 3);

    double xdisp = new_com[0] - old_com[0];
    double ydisp = new_com[1] - old_com[1];
    double zdisp = new_com[2] - old_com[2];
    
    // std::array<double, 9> new_quadrupole{};
    // for (int j = 0; j < 3; j++) {
    //   for (int i = 0; i < 3; i++) {
    //     new_quadrupole[3*i + j] = old_quadrupole[3*i + j] + tot_mass * disp[i] * disp[j];
    //   }
    // }

    std::array<double, 6> new_quadrupole{};
    // for (int j = 0; j < 3; j++) {
    //   for (int i = 0; i <= j; i++) {
    //     new_quadrupole[3*i + j] = old_quadrupole[3*i + j] + tot_mass * disp[i] * disp[j];
    //   }
    // }

    new_quadrupole[0] = old_quadrupole[0] + tot_mass * xdisp * xdisp;
    new_quadrupole[1] = old_quadrupole[1] + tot_mass * xdisp * ydisp;
    new_quadrupole[2] = old_quadrupole[2] + tot_mass * xdisp * zdisp;
    new_quadrupole[3] = old_quadrupole[3] + tot_mass * ydisp * ydisp;
    new_quadrupole[4] = old_quadrupole[4] + tot_mass * ydisp * zdisp;
    new_quadrupole[5] = old_quadrupole[5] + tot_mass * zdisp * zdisp;
    
    return new_quadrupole;

  }

  // subtract two arrays
  // we only ever use this with vectors (size=3)
  // any way to do this with arrays? size is not known at compile time
  // will probably just get rid of this
  static std::vector<double> subtract_(std::vector<double> a, std::vector<double> b, int size) throw()
  {
    // should assert that a and b have size "size"
    std::vector<double> diff (size, 0);
    
    for (int i = 0; i < size; i++) {
      diff[i] = a[i] - b[i];
    }

    return diff;
    
  }

  // add two arrays
  // will need to create add_11, add_22, add_33
  static std::vector<double> add_(std::vector<double> a, std::vector<double> b, int size) throw()
  {
    // should assert that a and b have size "size"
    std::vector<double> sum (size, 0);

    for (int i = 0; i < size; i++) {
      sum[i] = a[i] + b[i];
    }

    return sum;
  }

  // compute the outer product of two (3-elements) vectors 
  static std::vector<double> outer_11_(std::vector<double> a, std::vector<double> b) throw()
  {
    // assert that a and b are 3-element vectors 
    std::vector<double> prod (6, 0);

    // for (int j = 0; j < 3; j++) {
    //   for (int i = 0; i < 3; i++) {
    //     prod[3*i + j] = a[i] * b[j];
    //   }
    // }

    // we only ever have to take the outer product of a vector with itself, so should change arguments
    // for (int j = 0; j < 3; j++) {
    //   for (int i = 0; i <= j; i++) {
    //     prod[3*i + j] = a[i] * b[j];
    //   }
    // }

    prod[0] = a[0] * b[0];
    prod[1] = a[0] * b[1];
    prod[2] = a[0] * b[2];
    prod[3] = a[1] * b[1];
    prod[4] = a[1] * b[2];
    prod[5] = a[2] * b[2];

    return prod;
  }

  // compute the outer product of a vector 'a' with a rank-2 tensor 'b'
  static std::vector<double> outer_12_(std::vector<double> a, std::vector<double> b) throw()
  {
    // std::vector<double> prod (27, 0);

    // for (int k = 0; k < 3; k++) {
    //   for (int j = 0; j < 3; j++) {
    //     for (int i = 0; i < 3; i++) {
    //       prod[9*k + 3*i + j] = a[i] * b[3*j + k];
    //     }
    //   }
    // }

    // as with outer_11, outer_12 is only ever symmetric in our cases
    std::vector<double> prod (10, 0);

    // for (int k = 0; k < 3; k++) {
    //   for (int j = 0; j <= k; j++) {
    //     for (int i = 0; i <= j; i++) {
    //       prod[9*k + 3*i + j] = a[i] * b[3*j + k];
    //     }
    //   }
    // }

    prod[0] = a[0] * b[0];  // 000
    prod[1] = a[0] * b[1];  // 001
    prod[2] = a[0] * b[2];  // 002
    prod[3] = a[0] * b[3];  // 011
    prod[4] = a[0] * b[4];  // 012
    prod[5] = a[0] * b[5];  // 022
    prod[6] = a[1] * b[3];  // 111
    prod[7] = a[1] * b[4];  // 112
    prod[8] = a[1] * b[5];  // 122
    prod[9] = a[2] * b[5];  // 222


    return prod;
  }

  // multiply an array by a scalar
  // will need to create multiple versions of this for sizes 3, 6, 10
  static std::vector<double> dot_scalar_(double a, std::vector<double> b, int size) throw()
  {
    // double * prod = new double[size];
    std::vector<double> prod (size, 0);

    for (int i = 0; i < size; i++) {
      prod[i] = a * b[i];
    }

    return prod;
  }

  // compute the dot product of two vectors
  static double dot_11_(std::vector<double> a, std::vector<double> b) throw()
  {
    // double prod = 0;

    // for (int i = 0; i < 3; i++) {
    //   prod += a[i] * b[i];
    // }

    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  // compute the dot product of a vector 'a' with a rank-2 tensor 'b'
  static std::vector<double> dot_12_(std::vector<double> a, std::vector<double> b) throw()
  {
    // double * prod = new double[3];
    std::vector<double> prod (3, 0);

    // for (int j = 0; j < 3; j++) {
    //   for (int i = 0; i < 3; i++) {
    //     prod[j] += a[i] * b[3*i + j];
    //   }
    // }

    prod[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];   // 00 + 01 + 02
    prod[1] = a[0]*b[1] + a[1]*b[3] + a[2]*b[4];   // 01 + 11 + 12
    prod[2] = a[0]*b[2] + a[1]*b[4] + a[2]*b[5];   // 02 + 12 + 22
    

    return prod;
  }

  // compute the dot product of a vector 'a' with a rank-3 tensor 'b'
  static std::vector<double> dot_13_(std::vector<double> a, std::vector<double> b) throw()
  {

    // double * prod = new double[9];
    std::vector<double> prod (6, 0);

    // for (int k = 0; k < 3; k++) {
    //   for (int j = 0; j < 3; j++) {
    //     for (int i = 0; i < 3; i++) {
    //       prod[3*j + k] += a[i] * b[9*k + 3*i + j];
    //     }
    //   }
    // }

    prod[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    prod[1] = a[0]*b[1] + a[1]*b[3] + a[2]*b[4];
    prod[2] = a[0]*b[2] + a[1]*b[4] + a[2]*b[5];
    prod[3] = a[0]*b[3] + a[1]*b[6] + a[2]*b[7];
    prod[4] = a[0]*b[4] + a[1]*b[7] + a[2]*b[8];
    prod[5] = a[0]*b[5] + a[1]*b[8] + a[2]*b[9];

    return prod;

  }

  // compute the dot product of a rank-1 tensor 'a' with a rank-4 tensor 'b'
  static std::vector<double> dot_14_(std::vector<double> a, std::vector<double> b) throw()
  {
    std::vector<double> prod (10, 0);

    // for (int l = 0; l < 3; l++) {
    //   for (int k = 0; k < 3; k++) {
    //     for (int j = 0; j < 3; j++) {
    //       for (int i = 0; i < 3; i++) {
    //         prod[k + 3*(j + 3*l)] += a[i] * b[27*l + 9*k + 3*i + j];
    //       }
    //     }
    //   }
    // }

    prod[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    prod[1] = a[0]*b[1] + a[1]*b[3] + a[2]*b[4];
    prod[2] = a[0]*b[2] + a[1]*b[4] + a[2]*b[5];
    prod[3] = a[0]*b[3] + a[1]*b[6] + a[2]*b[7];
    prod[4] = a[0]*b[4] + a[1]*b[7] + a[2]*b[8];
    prod[5] = a[0]*b[5] + a[1]*b[8] + a[2]*b[9];
    prod[6] = a[0]*b[6] + a[1]*b[10] + a[2]*b[11];
    prod[7] = a[0]*b[7] + a[1]*b[11] + a[2]*b[12];
    prod[8] = a[0]*b[8] + a[1]*b[12] + a[2]*b[13];
    prod[9] = a[0]*b[9] + a[1]*b[13] + a[2]*b[14];

    return prod;
  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-2 tensor 'b'
  // don't really need this (only for d0)
  static double dot_22_(std::vector<double> a, std::vector<double> b) throw()
  {
    // double prod = 0;
    
    // for (int j = 0; j < 3; j++) {
    //   for (int i = 0; i < 3; i++) {
    //     prod += a[3*i + j] * b[3*i + j];
    //   }
    // }

    return (a[0]*b[0] + 2*a[1]*b[1] + 2*a[2]*b[2] 
           + a[3]*b[3] + 2*a[4]*b[4] 
           + a[5]*b[5]);
  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-3 tensor 'b'
  static std::vector<double> dot_23_(std::vector<double> a, std::vector<double> b) throw()
  {
    std::vector<double> prod (3, 0);

    // for (int k = 0; k < 3; k++) {
    //   for (int j = 0; j < 3; j++) {
    //     for (int i = 0; i < 3; i++) {
    //       prod[k] += a[3*i + j] * b[9*k + 3*i + j];
    //     }
    //   }
    // }

    prod[0] = a[0]*b[0] + 2*a[1]*b[1] + 2*a[2]*b[2] + a[3]*b[3] + 2*a[4]*b[4] + a[5]*b[5];
    prod[1] = a[0]*b[1] + 2*a[1]*b[3] + 2*a[2]*b[4] + a[3]*b[6] + 2*a[4]*b[7] + a[5]*b[8];
    prod[2] = a[0]*b[2] + 2*a[1]*b[4] + 2*a[2]*b[5] + a[3]*b[7] + 2*a[4]*b[8] + a[5]*b[9];
 
    return prod;
  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-4 tensor 'b'
  static std::vector<double> dot_24_(std::vector<double> a, std::vector<double> b) throw()
  {
    std::vector<double> prod (6, 0);

    // for (int l = 0; l < 3; l++) {
    //   for (int k = 0; k < 3; k++) {
    //     for (int j = 0; j < 3; j++) {
    //       for (int i = 0; i < 3; i++) {
    //         prod[j + 3*i] += a[3*i + j] * b[27*l + 9*k + 3*i + j];
    //       }
    //     }
    //   }
    // }

    prod[0] = a[0]*b[0] + 2*a[1]*b[1] + 2*a[2]*b[2] + a[3]*b[3] + 2*a[4]*b[4] + a[5]*b[5];
    prod[1] = a[0]*b[1] + 2*a[1]*b[3] + 2*a[2]*b[4] + a[3]*b[6] + 2*a[4]*b[7] + a[5]*b[8];
    prod[2] = a[0]*b[2] + 2*a[1]*b[4] + 2*a[2]*b[5] + a[3]*b[7] + 2*a[4]*b[8] + a[5]*b[9];
    prod[3] = a[0]*b[3] + 2*a[1]*b[6] + 2*a[2]*b[7] + a[3]*b[10] + 2*a[4]*b[11] + a[5]*b[12];
    prod[4] = a[0]*b[4] + 2*a[1]*b[7] + 2*a[2]*b[8] + a[3]*b[11] + 2*a[4]*b[12] + a[5]*b[13];
    prod[5] = a[0]*b[5] + 2*a[1]*b[8] + 2*a[2]*b[9] + a[3]*b[12] + 2*a[4]*b[13] + a[5]*b[14];

    return prod;
  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-5 tensor 'b'
  static std::vector<double> dot_25_(std::vector<double> a, std::vector<double> b) throw()
  {
    std::vector<double> prod (10, 0);

    // for (int m = 0; m < 3; m++) {
    //   for (int l = 0; l < 3; l++) {
    //     for (int k = 0; k < 3; k++) {
    //       for (int j = 0; j < 3; j++) {
    //         for (int i = 0; i < 3; i++) {
    //           prod[9*k + 3*i + j] += a[3*i + j] * b[81*m + 27*l + 9*k + 3*i + j];
    //         }
    //       }
    //     }
    //   }
    // }

    prod[0] = a[0]*b[0] + 2*a[1]*b[1] + 2*a[2]*b[2] + a[3]*b[3] + 2*a[4]*b[4] + a[5]*b[5];
    prod[1] = a[0]*b[1] + 2*a[1]*b[3] + 2*a[2]*b[4] + a[3]*b[6] + 2*a[4]*b[7] + a[5]*b[8];
    prod[2] = a[0]*b[2] + 2*a[1]*b[4] + 2*a[2]*b[5] + a[3]*b[7] + 2*a[4]*b[8] + a[5]*b[9];
    prod[3] = a[0]*b[3] + 2*a[1]*b[6] + 2*a[2]*b[7] + a[3]*b[10] + 2*a[4]*b[11] + a[5]*b[12];
    prod[4] = a[0]*b[4] + 2*a[1]*b[7] + 2*a[2]*b[8] + a[3]*b[11] + 2*a[4]*b[12] + a[5]*b[13];
    prod[5] = a[0]*b[5] + 2*a[1]*b[8] + 2*a[2]*b[9] + a[3]*b[12] + 2*a[4]*b[13] + a[5]*b[14];
    prod[6] = a[0]*b[6] + 2*a[1]*b[10] + 2*a[2]*b[11] + a[3]*b[15] + 2*a[4]*b[16] + a[5]*b[17];
    prod[7] = a[0]*b[7] + 2*a[1]*b[11] + 2*a[2]*b[12] + a[3]*b[16] + 2*a[4]*b[17] + a[5]*b[18];
    prod[8] = a[0]*b[8] + 2*a[1]*b[12] + 2*a[2]*b[13] + a[3]*b[17] + 2*a[4]*b[18] + a[5]*b[19];
    prod[9] = a[0]*b[9] + 2*a[1]*b[13] + 2*a[2]*b[14] + a[3]*b[18] + 2*a[4]*b[19] + a[5]*b[20];

    return prod;
  }

  // compute the dot product of a rank-3 tensor 'a' with a rank-3 tensor 'b'
  // don't really need this (only for d0)
  static double dot_33_(std::vector<double> a, std::vector<double> b) throw()
  {
    // double prod = 0;

    // for (int k = 0; k < 3; k++) {
    //   for (int j = 0; j < 3; j++) {
    //     for (int i = 0; i < 3; i++) {
    //       prod += a[9*k + 3*i + j] * b[9*k + 3*i + j];
    //     }
    //   }
    // }

    return (a[0]*b[0] + 3*a[1]*b[1] + 3*a[2]*b[2]
           + 3*a[3]*b[3] + 6*a[4]*b[4] + 3*a[5]*b[5]
           + a[6]*b[6] + 3*a[7]*b[7] + 3*a[8]*b[8] + a[9]*b[9]);
  }

  // compute the dot product of a rank-3 tensor 'a' with a rank-4 tensor 'b'
  static std::vector<double> dot_34_(std::vector<double> a, std::vector<double> b) throw()
  {
    std::vector<double> prod (3, 0);

    // for (int l = 0; l < 3; l++) {
    //   for (int k = 0; k < 3; k++) {
    //     for (int j = 0; j < 3; j++) {
    //       for (int i = 0; i < 3; i++) {
    //         prod[l] += a[9*k + 3*i + j] * b[27*l + 9*k + 3*i + j];
    //       }
    //     }
    //   }
    // }

    prod[0] = a[0]*b[0] + 3*a[1]*b[1] + 3*a[2]*b[2] + 3*a[3]*b[3] + 6*a[4]*b[4] + 3*a[5]*b[5] + a[6]*b[6] + 3*a[7]*b[7] + 3*a[8]*b[8] + a[9]*b[9];
    prod[1] = a[0]*b[1] + 3*a[1]*b[3] + 3*a[2]*b[4] + 3*a[3]*b[6] + 6*a[4]*b[7] + 3*a[5]*b[8] + a[6]*b[10] + 3*a[7]*b[11] + 3*a[8]*b[12] + a[9]*b[13];
    prod[2] = a[0]*b[2] + 3*a[1]*b[4] + 3*a[2]*b[5] + 3*a[3]*b[7] + 6*a[4]*b[8] + 3*a[5]*b[9] + a[6]*b[11] + 3*a[7]*b[12] + 3*a[8]*b[13] + a[9]*b[14];

    return prod;
  }

  // compute the dot product of a rank-3 tensor 'a' with a rank-5 tensor 'b'
  static std::vector<double> dot_35_(std::vector<double> a, std::vector<double> b) throw()
  {
    std::vector<double> prod (6, 0);

    // for (int m = 0; m < 3; m++) {
    //   for (int l = 0; l < 3; l++) {
    //     for (int k = 0; k < 3; k++) {
    //       for (int j = 0; j < 3; j++) {
    //         for (int i = 0; i < 3; i++) {
    //           prod[3*l + m] += a[9*k + 3*i + j] * b[81*m + 27*l + 9*k + 3*i + j];
    //         }
    //       }
    //     }
    //   }
    // }

    prod[0] = a[0]*b[0] + 3*a[1]*b[1] + 3*a[2]*b[2] + 3*a[3]*b[3] + 6*a[4]*b[4] + 3*a[5]*b[5] + a[6]*b[6] + 3*a[7]*b[7] + 3*a[8]*b[8] + a[9]*b[9];
    prod[1] = a[0]*b[1] + 3*a[1]*b[3] + 3*a[2]*b[4] + 3*a[3]*b[6] + 6*a[4]*b[7] + 3*a[5]*b[8] + a[6]*b[10] + 3*a[7]*b[11] + 3*a[8]*b[12] + a[9]*b[13];
    prod[2] = a[0]*b[2] + 3*a[1]*b[4] + 3*a[2]*b[5] + 3*a[3]*b[7] + 6*a[4]*b[8] + 3*a[5]*b[9] + a[6]*b[11] + 3*a[7]*b[12] + 3*a[8]*b[13] + a[9]*b[14];
    prod[3] = a[0]*b[3] + 3*a[1]*b[6] + 3*a[2]*b[7] + 3*a[3]*b[10] + 6*a[4]*b[11] + 3*a[5]*b[12] + a[6]*b[15] + 3*a[7]*b[16] + 3*a[8]*b[17] + a[9]*b[18];
    prod[4] = a[0]*b[4] + 3*a[1]*b[7] + 3*a[2]*b[8] + 3*a[3]*b[11] + 6*a[4]*b[12] + 3*a[5]*b[13] + a[6]*b[16] + 3*a[7]*b[17] + 3*a[8]*b[18] + a[9]*b[19];
    prod[5] = a[0]*b[5] + 3*a[1]*b[8] + 3*a[2]*b[9] + 3*a[3]*b[12] + 6*a[4]*b[13] + 3*a[5]*b[14] + a[6]*b[17] + 3*a[7]*b[18] + 3*a[8]*b[19] + a[9]*b[20];


    return prod;
  }

  // compute the dot product of a rank-3 tensor 'a' with a rank-6 tensor 'b'
  static std::vector<double> dot_36_(std::vector<double> a, std::vector<double> b) throw()
  {
    std::vector<double> prod (10, 0);
    
    // for (int n = 0; n < 3; n++) {
    //   for (int m = 0; m < 3; m++) {
    //     for (int l = 0; l < 3; l++) {
    //       for (int k = 0; k < 3; k++) {
    //         for (int j = 0; j < 3; j++) {
    //           for (int i = 0; i < 3; i++) {
    //             prod[9*n + 3*l + m] += a[9*k + 3*i + j] * b[243*n + 81*m + 27*l + 9*k + 3*i + j];
    //           }
    //         }
    //       }
    //     }
    //   }
    // }

    prod[0] = a[0]*b[0] + 3*a[1]*b[1] + 3*a[2]*b[2] + 3*a[3]*b[3] + 6*a[4]*b[4] + 3*a[5]*b[5] + a[6]*b[6] + 3*a[7]*b[7] + 3*a[8]*b[8] + a[9]*b[9];
    prod[1] = a[0]*b[1] + 3*a[1]*b[3] + 3*a[2]*b[4] + 3*a[3]*b[6] + 6*a[4]*b[7] + 3*a[5]*b[8] + a[6]*b[10] + 3*a[7]*b[11] + 3*a[8]*b[12] + a[9]*b[13];
    prod[2] = a[0]*b[2] + 3*a[1]*b[4] + 3*a[2]*b[5] + 3*a[3]*b[7] + 6*a[4]*b[8] + 3*a[5]*b[9] + a[6]*b[11] + 3*a[7]*b[12] + 3*a[8]*b[13] + a[9]*b[14];
    prod[3] = a[0]*b[3] + 3*a[1]*b[6] + 3*a[2]*b[7] + 3*a[3]*b[10] + 6*a[4]*b[11] + 3*a[5]*b[12] + a[6]*b[15] + 3*a[7]*b[16] + 3*a[8]*b[17] + a[9]*b[18];
    prod[4] = a[0]*b[4] + 3*a[1]*b[7] + 3*a[2]*b[8] + 3*a[3]*b[11] + 6*a[4]*b[12] + 3*a[5]*b[13] + a[6]*b[16] + 3*a[7]*b[17] + 3*a[8]*b[18] + a[9]*b[19];
    prod[5] = a[0]*b[5] + 3*a[1]*b[8] + 3*a[2]*b[9] + 3*a[3]*b[12] + 6*a[4]*b[13] + 3*a[5]*b[14] + a[6]*b[17] + 3*a[7]*b[18] + 3*a[8]*b[19] + a[9]*b[20];
    prod[6] = a[0]*b[6] + 3*a[1]*b[10] + 3*a[2]*b[11] + 3*a[3]*b[15] + 6*a[4]*b[16] + 3*a[5]*b[17] + a[6]*b[21] + 3*a[7]*b[22] + 3*a[8]*b[23] + a[9]*b[24];
    prod[7] = a[0]*b[7] + 3*a[1]*b[11] + 3*a[2]*b[12] + 3*a[3]*b[16] + 6*a[4]*b[17] + 3*a[5]*b[18] + a[6]*b[22] + 3*a[7]*b[23] + 3*a[8]*b[24] + a[9]*b[25];
    prod[8] = a[0]*b[8] + 3*a[1]*b[12] + 3*a[2]*b[13] + 3*a[3]*b[17] + 6*a[4]*b[18] + 3*a[5]*b[19] + a[6]*b[23] + 3*a[7]*b[24] + 3*a[8]*b[25] + a[9]*b[26];
    prod[9] = a[0]*b[9] + 3*a[1]*b[13] + 3*a[2]*b[14] + 3*a[3]*b[18] + 6*a[4]*b[19] + 3*a[5]*b[20] + a[6]*b[24] + 3*a[7]*b[25] + 3*a[8]*b[26] + a[9]*b[27];


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

  /// Gravitational softening length
  double eps0_;

  /// Cutoff distance for gravitational softening
  double r0_;

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

  /// dimensions of interpolation grid for periodic Ewald summation
  int interp_xpoints_;
  int interp_ypoints_;
  int interp_zpoints_;

  EnzoMethodEwald * ewald_;

};

#endif /* ENZO_ENZO_METHOD_MULTIPOLE_HPP */
