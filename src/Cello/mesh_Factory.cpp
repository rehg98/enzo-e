// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#include "mesh.hpp"
#include "charm_simulation.hpp"

//----------------------------------------------------------------------

Hierarchy * Factory::create_hierarchy
( int refinement, int min_level, int max_level) const throw ()
{
  return new Hierarchy (this,refinement, min_level, max_level);
}

//----------------------------------------------------------------------

void Factory::pup (PUP::er &p)

{
  TRACEPUP;

  // NOTE: change this function whenever attributes change

  PUP::able::pup(p);

  //  p | bound_arrays_;

}

//----------------------------------------------------------------------

IoBlock * Factory::create_io_block () const throw()
{
  return new IoBlock;
}

//----------------------------------------------------------------------

IoFieldData * Factory::create_io_field_data () const throw()
{
  return new IoFieldData();

}

//----------------------------------------------------------------------

IoParticleData * Factory::create_io_particle_data () const throw()
{
  return new IoParticleData();
}

//----------------------------------------------------------------------

CProxy_Block Factory::new_block_proxy
(
 DataMsg * data_msg,
 int nbx, int nby, int nbz
 ) const throw()
{
  TRACE3("Factory::new_block_array(na(%d %d %d)))",
	  nbx,nby,nbz);

  CProxy_Block proxy_block;

  CProxy_MappingArray array_map  = CProxy_MappingArray::ckNew(nbx,nby,nbz);

  CkArrayOptions opts;
  opts.setMap(array_map);
  proxy_block = CProxy_Block::ckNew(opts);

  return proxy_block;
}

//----------------------------------------------------------------------

void Factory::create_block_array
(
 DataMsg * data_msg,
 CProxy_Block proxy_block,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_data
 ) const throw()
{
  TRACE7("Factory::create_block_array(na(%d %d %d) n(%d %d %d) num_field_data %d)",
	 nbx,nby,nbz,nx,ny,nz,num_field_data);

  int count_adapt;

  int    cycle = 0;
  double time  = 0.0;
  double dt    = 0.0;
  int num_face_level = 0;
  int * face_level = 0;

  for (int ix=0; ix<nbx; ix++) {
    for (int iy=0; iy<nby; iy++) {
      for (int iz=0; iz<nbz; iz++) {

	Index index(ix,iy,iz);

	MsgRefine * msg = new MsgRefine
	  (index,
	   nx,ny,nz,
	   num_field_data,
	   count_adapt = 0,
	   cycle,time,dt,
	   refresh_same,
	   num_face_level, face_level, nullptr);

	msg->set_data_msg(data_msg);
#ifdef BUG_FIX_150
	cello::simulation()->set_msg_refine (index,msg);
	proxy_block[index].insert (process_type(CkMyPe()));
#else
	proxy_block[index].insert (msg);
#endif
	// --------------------------------------------------

      }
    }
  }

  TRACE1("Factory::create_block_array = %p",&proxy_block);
}

//----------------------------------------------------------------------

void Factory::create_subblock_array
(
 DataMsg * data_msg,
 CProxy_Block block_array,
 int min_level,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks
 ) const throw()
{
  TRACE8("Factory::create_subblock_array(min_level %d na(%d %d %d) n(%d %d %d) num_field_blocks %d",
	 min_level,nbx,nby,nbz,nx,ny,nz,num_field_blocks);

  if (min_level >= 0) {
    WARNING1("Factor::create_subblock_array",
	     "Trying to create subblock array with min_level %d >= 0",
	     min_level);
    return ;
  }

  for (int level = -1; level >= min_level; level--) {

    if (nbx > 1) nbx = ceil(0.5*nbx);
    if (nby > 1) nby = ceil(0.5*nby);
    if (nbz > 1) nbz = ceil(0.5*nbz);

    int count_adapt;

    int    cycle = 0;
    double time  = 0.0;
    double dt    = 0.0;
    int num_face_level = 0;
    int * face_level = 0;

    for (int ix=0; ix<nbx; ix++) {
      for (int iy=0; iy<nby; iy++) {
	for (int iz=0; iz<nbz; iz++) {

	  int shift = -level;

	  Index index(ix<<shift,iy<<shift,iz<<shift);

	  index.set_level(level);

	  TRACE3 ("inserting %d %d %d",ix,iy,iz);

	  MsgRefine * msg = new MsgRefine
	    (index,
	     nx,ny,nz,
	     num_field_blocks,
	     count_adapt = 0,
	     cycle,time,dt,
	     refresh_same,
	     num_face_level, face_level, nullptr);

	  msg->set_data_msg(data_msg);

#ifdef BUG_FIX_150
          cello::simulation()->set_msg_refine (index,msg);
          block_array[index].insert (process_type(CkMyPe()));
#else
          block_array[index].insert (msg);
#endif
	  // --------------------------------------------------

	}
      }
    }
  }
}

//----------------------------------------------------------------------

void Factory::create_block
(
 DataMsg * data_msg,
 CProxy_Block block_array,
 Index index,
 int nx, int ny, int nz,
 int num_field_data,
 int count_adapt,
 int cycle, double time, double dt,
 int narray, char * array, int refresh_type,
#ifdef OLD_ADAPT
 int num_face_level,
 int * face_level,
#endif
 Adapt * adapt,
 Simulation * simulation
 ) const throw()
{

  TRACE3("Factory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("Factory::create_block(num_field_data %d  count_adapt %d)",
	 num_field_data,count_adapt);

  // --------------------------------------------------
  // ENTRY: #3 Factory::create_block() -> Block::Block()
  // ENTRY: level > 0 block array insert
  // --------------------------------------------------

#ifdef NEW_ADAPT
  const int num_face_level = 27;
  int * face_level = new int[num_face_level];
  std::fill_n(face_level,num_face_level,INDEX_UNDEFINED_LEVEL);
  // for (int i=1; i<adapt->num_neighbors(); i++) {
  //   Index index_neighbor = adapt->index(i);
  // }
  // int * face_level = new int[27];
  // int level = index.level();
  // int cx,cy,cz;
  // index.child(level,&cx,&cy,&cz);
  // for (int iz=0; iz<3; iz++) {
  //   int ipz=(iz+1+cz) / 2;
  //   for (int iy=0; iy<3; iy++) {
  //     int ipy=(iy+1+cy) / 2;
  //     for (int ix=0; ix<3; ix++) {
  //       int ipx=(ix+1+cx) / 2;
  //       int i=ix+3*(iy+3*iz);
  //       int ip=ipx+3*(ipy+3*ipz);
  //       face_level[i] = face_level_parent[ip];
  //     }
  //   }
  // }
#endif

  MsgRefine * msg = new MsgRefine
    (index,
     nx,ny,nz,
     num_field_data,
     count_adapt,
     cycle,time,dt,
     refresh_type,
     num_face_level, face_level, adapt
     );

  msg->set_data_msg (data_msg);

#ifdef BUG_FIX_150
  cello::simulation()->set_msg_refine (index,msg);
  block_array[index].insert (process_type(CkMyPe()));
#else
  block_array[index].insert (msg);
#endif
}

