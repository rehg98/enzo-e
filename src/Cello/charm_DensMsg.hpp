// See LICENSE_CELLO file for license and copyright information

/// @file     charm_DensMsg.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     2022-11-22
/// @brief    [\ref Charm] Declaration of the DensMsg Charm++ Message

#ifndef CHARM_DENS_MSG_HPP
#define CHARM_DENS_MSG_HPP

class DensMsg : public CMessage_DensMsg {

public: // attributes

  /// lower left corner of Block
  double lo[3];

  std::vector<double> dens;

};

#endif /* CHARM_DENS_MSG_HPP */
