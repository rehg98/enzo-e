// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMethodEwald.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu) 
/// @date     Thursday September 28, 2023
/// @brief    [\ref Enzo] Declaration of EnzoMethodEwald
///           Compute Ewald sums for periodic boundary conditions

#ifndef ENZO_METHOD_EWALD_HPP
#define ENZO_METHOD_EWALD_HPP

#include "view.hpp"

class EnzoMethodEwald {

  /// @class    EnzoMethodEwald
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Compute Ewald sums for periodic boundary conditions

public: 

  /// Create a new EnzoMethodEwald object
  EnzoMethodEwald (int interp_xpoints, int interp_ypoints, int interp_zpoints);

  EnzoMethodEwald()
    : interp_xpoints_(64),  // number of interpolation points in the x-direction
      interp_ypoints_(64),  // number of interpolation points in the y-direction
      interp_zpoints_(64)   // number of interpolation points in the z-direction
  { }

  void pup(PUP::er &p) {
    p | interp_xpoints_;
    p | interp_ypoints_;
    p | interp_zpoints_;
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

  // double interp_d0(double x, double y, double z) throw();  // this function is not necessary for computing acceleration
  std::array<double, 3> interp_d1(double x, double y, double z) throw();
  std::array<double, 6> interp_d2(double x, double y, double z) throw();
  std::array<double, 10> interp_d3(double x, double y, double z) throw();

  // double d0(double x, double y, double z) throw();  // this function is not necessary for computing acceleration 
  // std::array<double, 3> d1(double x, double y, double z) throw();
  // std::array<double, 6> d2(double x, double y, double z) throw();
  // std::array<double, 10> d3(double x, double y, double z) throw();
  // std::array<double, 15> d4(double x, double y, double z) throw();
  // std::array<double, 21> d5(double x, double y, double z) throw();
  // std::array<double, 28> d6(double x, double y, double z) throw();

  CelloView<double, 1> d1(double x, double y, double z) throw();
  CelloView<double, 1> d2(double x, double y, double z) throw();
  CelloView<double, 1> d3(double x, double y, double z) throw();
  CelloView<double, 1> d4(double x, double y, double z) throw();
  CelloView<double, 1> d5(double x, double y, double z) throw();
  CelloView<double, 1> d6(double x, double y, double z) throw();

protected: //attributes

  // dimensions of interpolation grid (Nx x Ny x Nz)
  int interp_xpoints_;
  int interp_ypoints_;
  int interp_zpoints_;

  // int tot_points_ = interp_xpoints_ * interp_ypoints_ * interp_zpoints_;
  //const int tot_points_ = 10

  /// derivative interpolation arrays
  // std::array<double, tot_points_> d0_array_; 
  // std::array<std::array<double, 3>, tot_points_> d1_array_; 
  // std::array<std::array<double, 6>, tot_points_> d2_array_; 
  // std::array<std::array<double, 10>, tot_points_> d3_array_; 
  // std::array<std::array<double, 15>, tot_points_> d4_array_; 
  // std::array<std::array<double, 21>, tot_points_> d5_array_; 
  // std::array<std::array<double, 28>, tot_points_> d6_array_; 

  CelloView<double, 2> d1_array_; 
  CelloView<double, 2> d2_array_; 
  CelloView<double, 2> d3_array_; 
  CelloView<double, 2> d4_array_; 
  CelloView<double, 2> d5_array_; 
  CelloView<double, 2> d6_array_; 


  
  

};

#endif /* ENZO_METHOD_EWALD_HPP */
