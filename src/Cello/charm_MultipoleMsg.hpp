// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MultipoleMsg.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     2022-7-20
/// @brief    [\ref Charm] Declaration of the MultipoleMsg Charm++ Message

#ifndef CHARM_MULTIPOLE_MSG_HPP
#define CHARM_MULTIPOLE_MSG_HPP

class MultipoleMsg : public CMessage_MultipoleMsg {

public: // attributes

  int child_index() const { return ic3[0] + 2*(ic3[1] + 2*(ic3[2])); }
  
  /// total mass
  double mass;

  /// center of mass
  double com[3];

  /// quadrupole moment
  double quadrupole[9];

  /// Child indices
  int ic3[3];
};

#endif /* CHARM_MULTIPOLE_MSG_HPP */
