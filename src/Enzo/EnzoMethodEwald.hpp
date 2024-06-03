// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMethodEwald.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu) 
/// @date     Thursday September 28, 2023
/// @brief    [\ref Enzo] Declaration of EnzoMethodEwald
///           Compute Ewald sums for periodic boundary conditions

#ifndef ENZO_METHOD_EWALD_HPP
#define ENZO_METHOD_EWALD_HPP

class EnzoMethodEwald {

  /// @class    EnzoMethodEwald
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Compute Ewald sums for periodic boundary conditions

public: 

  /// Create a new EnzoMethodEwald object
  EnzoMethodEwald (int interp_xpoints, int interp_ypoints, int interp_zpoints);

  EnzoMethodEwald()
    : d0_array_(), // Nx x Ny x Nz x 1 (on down-sampled grid of dimension Nx x Ny x Nz)
      d1_array_(), // Nx x Ny x Nz x 3
      d2_array_(), // Nx x Ny x Nz x 9
      d3_array_(), // Nx x Ny x Nz x 27
      d4_array_(), // Nx x Ny x Nz x 81
      d5_array_(), // Nx x Ny x Nz x 243
      d6_array_(), // Nx x Ny x Nz x 729
      interp_xpoints_(64),  // number of interpolation points in the x-direction
      interp_ypoints_(64),  // number of interpolation points in the y-direction
      interp_zpoints_(64)   // number of interpolation points in the z-direction
      // is_init_(0)
  { }

  void pup(PUP::er &p) {
    p | d0_array_;
    p | d1_array_;
    p | d2_array_;
    p | d3_array_;
    p | d4_array_;
    p | d5_array_;
    p | d6_array_;
    p | interp_xpoints_;
    p | interp_ypoints_;
    p | interp_zpoints_;
    // p | is_init_;
  }

  void init_interpolate_() throw();

  void find_nearest_interp_point(double x, double y, double z, 
                                 double * interp_x, double * interp_y, double * interp_z, int * i) throw()
  {
    double lox, loy, loz; 
    double hix, hiy, hiz;
    cello::hierarchy()->lower(&lox, &loy, &loz);
    cello::hierarchy()->upper(&hix, &hiy, &hiz);

    double dx = (hix - lox) / (interp_xpoints_ - 1);
    double dy = (hiy - loy) / (interp_ypoints_ - 1);
    double dz = (hiz - loz) / (interp_zpoints_ - 1);

    int ix = (int)((x-lox)/dx);
    int iy = (int)((y-loy)/dy);
    int iz = (int)((z-loz)/dz);

    *i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);

    *interp_x = lox + ix*dx;
    *interp_y = loy + iy*dy;
    *interp_z = loz + iz*dz;

    CkPrintf("%f, %f, %f, %d, %f, %f, %f\n", x, y, z, *i, *interp_x, *interp_y, *interp_z);

    return;
  }

  double interp_d0(double x, double y, double z) throw();  // this function is not necessary
  std::vector<double> interp_d1(double x, double y, double z) throw();
  std::vector<double> interp_d2(double x, double y, double z) throw();
  std::vector<double> interp_d3(double x, double y, double z) throw();

  double d0(double x, double y, double z) throw();  // this function is not necessary
  std::vector<double> d1(double x, double y, double z) throw();
  std::vector<double> d2(double x, double y, double z) throw();
  std::vector<double> d3(double x, double y, double z) throw();
  std::vector<double> d4(double x, double y, double z) throw();
  std::vector<double> d5(double x, double y, double z) throw();
  std::vector<double> d6(double x, double y, double z) throw();


protected: // attributes

  /// derivative interpolation arrays
  std::vector<double> d0_array_; 
  std::vector<std::vector<double>> d1_array_; 
  std::vector<std::vector<double>> d2_array_; 
  std::vector<std::vector<double>> d3_array_; 
  std::vector<std::vector<double>> d4_array_; 
  std::vector<std::vector<double>> d5_array_; 
  std::vector<std::vector<double>> d6_array_; 
  
  // dimensions of interpolation grid (Nx x Ny x Nz)
  int interp_xpoints_;
  int interp_ypoints_;
  int interp_zpoints_;

  // bool is_init_; // is_initialized ==> change from 0 to false

};

#endif /* ENZO_METHOD_EWALD_HPP */
