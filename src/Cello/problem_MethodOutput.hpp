// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOutput.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2021-03-09
/// @brief    [\ref Problem] Declaration for the MethodOutput class

#ifndef PROBLEM_METHOD_OUTPUT_HPP
#define PROBLEM_METHOD_OUTPUT_HPP

class Io;
class IoBlock;
class IoFieldData;
class IoParticleData;
class FileHdf5;

class MethodOutput : public Method
{
  /// @class    MethodOutput
  /// @ingroup  MethodOutput
  /// @brief    [\ref MethodOutput] Declaration of MethodOutput
  ///
  /// Method for writing data to disk files.  Designed to be highly
  /// scalable.

public: // interface

  /// Create a new MethodOutput
  MethodOutput
  (const Factory * factory,
   std::vector< std::string > file_name,
   std::vector< std::string > path_name,
   std::vector< std::string > field_list,
   std::vector< std::string > particle_list,
   int ghost_depth,
   int min_face_rank,
   bool all_fields,
   bool all_particles,
   int blocking_x = 1,
   int blocking_y = 1,
   int blocking_z = 1);

  /// Destructor
  virtual ~MethodOutput() throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodOutput);

  /// Charm++ PUP::able migration constructor
  MethodOutput (CkMigrateMessage *m)
    : Method(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void compute_continue (Block * block);
  void next (Block * block, MsgOutput *);
  void write (Block * block, MsgOutput *);

  const Factory * factory () const
  { return factory_; }

  public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute ( Block * block) throw();

  /// barrier before exiting
  void compute_done (Block * block);
  
  /// Return the name of this MethodOutput
  virtual std::string name () throw ()
  { return "output"; }

protected: // functions

  void output_ (Block * block);

  int is_writer_ (Index index);

  FileHdf5 * file_open_(Block * block, int a3[3]);
  void file_write_hierarchy_(FileHdf5 * file);
  void file_write_block_(FileHdf5 * , Block * , MsgOutput *);
  int file_count_(Block * block);
  void write_meta_ ( FileHdf5 * file, Io * io, std::string type_meta );
  DataMsg * create_data_msg_ (Block * block);
  
  
protected: // attributes

  /// File name and format
  std::vector <std::string> file_name_;

  /// Path name and format
  std::vector <std::string> path_name_;
  
  /// List of id's of fields to output
  std::vector<int> field_list_;

  /// List of id's of particle types to output
  std::vector<int> particle_list_;

  /// Ghost layer depth
  int ghost_depth_;

  /// Minimum dimensional face to output (0 corners, 1 edges, 2 facets)
  int min_face_rank_;

  /// Whether to output all fields, ignoring field_list_
  bool all_fields_;

  /// Whether to output all particles, ignoring particle_list_
  bool all_particles_;

  /// Size of the root-level octree array partitioning.  Data in all
  /// blocks in each partition are written to a single file by the
  /// minimal root-level Block
  int blocking_[3];

  /// Block Scalar int for counting files
  int is_count_;

  /// Factory for creating Io objects
  union {
    const Factory * factory_;
    // non-const USED IN pup() ONLY
    Factory * factory_nonconst_;
  };

};


#endif /* PROBLEM_METHOD_OUTPUT_HPP */
