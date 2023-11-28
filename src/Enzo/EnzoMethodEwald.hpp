// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMethodEwald.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu) 
/// @date     Thursday September 28, 2023
/// @brief    [\ref Enzo] Declaration of EnzoMethodEwald
///           Compute Ewald sums for periodic boundary conditions

#ifndef ENZO_METHOD_EWALD_HPP
#define ENZO_METHOD_EWALD_HPP

class EnzoMethodEwald : public Method {

  /// @class    EnzoMethodEwald
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Compute Ewald sums for periodic boundary conditions

public: 

  /// Create a new EnzoMethodEwald object
  EnzoMethodEwald (int interp_xpoints, int interp_ypoints, int interp_zpoints);

  EnzoMethodEwald()
    : Method(),
    d0_array_(), // Nx x Ny x Nz x 1 (on down-sampled grid of dimension Nx x Ny x Nz)
    d1_array_(), // Nx x Ny x Nz x 3
    d2_array_(), // Nx x Ny x Nz x 9
    d3_array_(), // Nx x Ny x Nz x 27
    d4_array_(), // Nx x Ny x Nz x 81
    d5_array_(), // Nx x Ny x Nz x 243
    d6_array_(), // Nx x Ny x Nz x 729
    interp_xpoints_(64),  // number of interpolation points in the x-direction
    interp_ypoints_(64),  // number of interpolation points in the y-direction
    interp_zpoints_(64)   // number of interpolation points in the z-direction

  { }


  void init_interpolate_() throw();

  void find_nearest_interp_point(double x, double y, double z, 
                                 double * interp_x, double * interp_y, double * interp_z, int * i) throw()
  {
    double lox, loy, loz; 
    double hix, hiy, hiz;
    cello::hierarchy()->lower(&lox, &loy, &loz);
    cello::hierarchy()->upper(&hix, &hiy, &hiz);

    double dx = (hix - lox) / interp_xpoints_;
    double dy = (hiy - loy) / interp_ypoints_;
    double dz = (hiz - loz) / interp_zpoints_;

    int ix = floor((x-lox)/dx);
    int iy = floor((y-loy)/dy);
    int iz = floor((z-loz)/dz);

    *i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);

    *interp_x = lox + (ix + 0.5)*dx;
    *interp_y = loy + (iy + 0.5)*dy;
    *interp_z = loz + (iz + 0.5)*dz;

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

};

#endif /* ENZO_METHOD_EWALD_HPP */