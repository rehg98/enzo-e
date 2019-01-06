#include "enzo.hpp"
#include "charm_enzo.hpp"


// This is a helper function for initializing the array from fields
void initialize_field_array(Block *block, EnzoArray<enzo_float> &array,
			    std::string field_name)
{
  Field field = block->data()->field();
  const int id = field.field_id(field_name);

  // get the field dimensions
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);
  enzo_float *data = (enzo_float *) field.values(field_name);
  array.initialize_wrapper(data, mz, my, mx);
}

// Class to be used for rotating axes.
// This class will be slightly extended and reused for implementing cosmic rays
// Using the coordinate systems given by eqn 68 of Gardiner & Stone
// (except that index numbers start at 0)
// The system of the mesh are x,y,z
// The system of the linear wave is x0, x1, x2
// The linear waves are defined along x0
// {{x0},   {{cos(a), 0, -sin(a)},   {{ cos(b), sin(b), 0},   {{x},
//  {x1}, =  {     0, 1,       0}, .  {-sin(b), cos(b), 0}, .  {y},
//  {x2}}    {sin(a), 0,  cos(a)}}    {      0,      0, 1}}    {z}}
//                    R2                          R1
//  R1, the second matix, rotates a 3D vector about it's third dimension by
//  angle -b. In this case, if z-axis points out of the page, then rotate
//  axes clockwise. The rotated y-axis now points along the x1 dimension
//
//  R2, the first matrix, rotates a 3D vector about it's second axis by angle
//  a. In this case, if the axis along the second dimension (x1) points out of
//  the page, then rotate the other axes counter clockwise.
//
//  Note that they give transformation of x0,x1,x2 -> x,y,z while the
//  calculations primarily use x,y,z -> x0,x1,x2 (since the the rotation is
//  easier to visualize

class Rotation {
public:
  Rotation(enzo_float a, enzo_float b)
  {
    matrix_[0][0] = std::cos(a)*std::cos(b);
    matrix_[0][1] = std::cos(a)*std::sin(b);
    matrix_[0][2] = std::sin(a);
    matrix_[1][0] = -1.*std::sin(b);
    matrix_[1][1] = std::cos(b);
    matrix_[1][2] = 0.;
    matrix_[2][0] = -1.*std::sin(a)*std::cos(b);
    matrix_[2][1] = -1.*std::sin(a)*std::sin(b);
    matrix_[2][2] = std::cos(a);
  }

  // rotates vector {v0, v1, v2} to {rot0, rot1, rot2}
  void rot(enzo_float v0, enzo_float v1, enzo_float v2,
	      enzo_float &rot0, enzo_float &rot1, enzo_float &rot2)
  {
    rot0 = matrix_[0][0]*v0 + matrix_[0][1]*v1 + matrix_[0][2];
    rot1 = matrix_[1][0]*v0 + matrix_[1][1]*v1 + matrix_[1][2];
    rot2 = matrix_[2][0]*v0 + matrix_[2][1]*v1 + matrix_[2][2];
  }

  // rotates vector {rot0, rot1, rot2} to {v0, v1, v2}
  void inv_rot(enzo_float rot0, enzo_float rot1, enzo_float rot2,
	       enzo_float &v0, enzo_float &v1, enzo_float &v2)
  {
    v0 = matrix_[0][0]*rot0 + matrix_[1][0]*rot1 + matrix_[2][0]*rot2;
    v1 = matrix_[0][1]*rot0 + matrix_[1][1]*rot1 + matrix_[2][1]*rot2;
    v2 = matrix_[0][2]*rot0 + matrix_[1][2]*rot1 + matrix_[2][2]*rot2;
  }

private:
  double matrix_[3][3];
};

// Use initializer functors to set up the grid
//   - for scalars: density, total energy density
//   - for vectors: momentum, magnetic vector potential
// These functors are most frequently most easily defined along axes that are
// rotated with respect to the grid. RotatedScalarInitializer and
// RotatedVectorInitializer are defined wrap around such initializers
// There are large similarities between Scalar and Vector Initializers, there
// may be a way to consolidate

class ScalarInitializer{
public:
  virtual ~ScalarInitializer()
  {}

  virtual enzo_float operator() (enzo_float x0, enzo_float x1, enzo_float x2)=0;
};

class LinearScalarInitializer : public ScalarInitializer{

public:
  LinearScalarInitializer(enzo_float background, enzo_float eigenvalue,
			  enzo_float amplitude, enzo_float lambda)
    : background_(background),
      eigenvalue_(eigenvalue),
      amplitude_(amplitude),
      lambda_(lambda)
  {}

  enzo_float operator() (enzo_float x0, enzo_float x1, enzo_float x2){
    return (background_ + amplitude_ * eigenvalue_
	    * std::cos(x0*2.*cello::pi/lambda_));
  }

protected:
  // value of the background state
  enzo_float background_;
  // value of the eigenvalue for a given mode
  enzo_float eigenvalue_;
  // perturbation amplitude
  enzo_float amplitude_;
  // wavelength
  enzo_float lambda_;
};


class RotatedScalarInitializer : public ScalarInitializer {
public:
  RotatedScalarInitializer(Rotation *rot, ScalarInitializer *inner)
    : ScalarInitializer(),
      rot_(rot),
      inner_(inner)
  {}

  ~RotatedScalarInitializer()
  {
    delete inner_;
    rot_ = NULL;
  }

  enzo_float operator() (enzo_float x, enzo_float y, enzo_float z)
  {
    enzo_float x0,x1,x2;
    rot_->rot(x,y,z,x0,x1,x2);
    return (*inner_)(x0,x1,x2);
  }

private:
  Rotation *rot_;
  ScalarInitializer *inner_;
};




class VectorInitializer {
public:
  virtual ~VectorInitializer()
  {}

  virtual void operator() (enzo_float x0, enzo_float x1, enzo_float x2,
			   enzo_float &v0, enzo_float &v1, enzo_float &v2)=0;
};

class LinearVectorInitializer : public VectorInitializer {

public:
  LinearVectorInitializer(enzo_float back0, enzo_float back1, enzo_float back2,
			  enzo_float ev0, enzo_float ev1, enzo_float ev2,
			  enzo_float amplitude, enzo_float lambda)
    : VectorInitializer(),
      back0_(back0), back1_(back1), back2_(back2),
      ev0_(ev0), ev1_(ev1), ev2_(ev2),
      amplitude_(amplitude),
      lambda_(lambda)
  {}

  void operator() (enzo_float x0, enzo_float x1, enzo_float x2,
		   enzo_float &v0, enzo_float &v1, enzo_float &v2)
  {
    v0 = (back0_ + amplitude_ * ev0_* std::cos(x0*2.*cello::pi/lambda_));
    v1 = (back1_ + amplitude_ * ev1_* std::cos(x1*2.*cello::pi/lambda_));
    v2 = (back2_ + amplitude_ * ev2_* std::cos(x2*2.*cello::pi/lambda_));
  }

protected:
  // value of the background state
  enzo_float back0_, back1_, back2_;
  // value of the eigenvalue for a given mode
  enzo_float ev0_, ev1_, ev2_;
  // perturbation amplitude
  enzo_float amplitude_;
  // wavelength
  enzo_float lambda_;
};


class RotatedVectorInitializer : public VectorInitializer {
public:
  RotatedVectorInitializer(Rotation *rot, VectorInitializer *inner)
    : VectorInitializer(),
      rot_(rot),
      inner_(inner)
  {}

  ~RotatedVectorInitializer()
  {
    delete inner_;
    rot_ = NULL;
  }

  void operator() (enzo_float x, enzo_float y, enzo_float z,
		   enzo_float &v0, enzo_float &v1, enzo_float &v2)
  {
    enzo_float x0,x1,x2, rot0, rot1, rot2;
    rot_->rot(x,y,z,x0,x1,x2);
    (*inner_)(x0,x1,x2,rot0,rot1,rot2);
    rot_->inv_rot(rot0,rot1,rot2,v0,v1,v2);
  }

private:
  Rotation *rot_;
  VectorInitializer *inner_;
};


class MeshPos {

public:
  MeshPos(EnzoBlock *block)
  {
    EnzoBlock * enzo_block = enzo::block(block);
    x_left_edge_ = block->GridLeftEdge[0];
    y_left_edge_ = block->GridLeftEdge[1];
    z_left_edge_ = block->GridLeftEdge[2];

    // ghost depth 
    gx_ = enzo_block->GridStartIndex[0];
    gy_ = enzo_block->GridStartIndex[1];
    gz_ = enzo_block->GridStartIndex[2];

    // cell widths
    dx_ = enzo_block->CellWidth[0];
    dy_ = enzo_block->CellWidth[1];
    dz_ = enzo_block->CellWidth[2];
  }

  // compute positions at cell-centers
  enzo_float x(int k, int j, int i){
    return x_left_edge_ + dx_* (0.5+(enzo_float)(i-gx_));}
  enzo_float y(int k, int j, int i){
    return y_left_edge_ + dy_* (0.5+(enzo_float)(j-gy_));}
  enzo_float z(int k, int j, int i){
    return z_left_edge_ + dz_* (0.5+(enzo_float)(k-gz_));}

  // computes the x, y, or z position for the cell corner at the
  // (i-1/2, j,k),  (i,j-1/2,k) or (i,j, k-1/2)
  enzo_float x_face(int k, int j, int i){
    return x_left_edge_ + dx_* (enzo_float)(i-gx_);}
  enzo_float y_face(int k, int j, int i){
    return y_left_edge_ + dy_* (enzo_float)(j-gy_);}
  enzo_float z_face(int k, int j, int i){
    return z_left_edge_ + dz_* (enzo_float)(k-gz_);}
  
private:
  // the starting position of the edge of the active grid
  enzo_float x_left_edge_, y_left_edge_, z_left_edge_;
  // ghost depth (and the index at the start of the active region)
  int gx_, gy_, gz_;
  enzo_float dx_, dy_, dz_;
};



void bfieldi_helper_(EnzoArray<enzo_float> &bfield, EnzoArray<enzo_float> &Aj,
		     EnzoArray<enzo_float> &Ak, int dim, enzo_float dj,
		     enzo_float dk)
{

  // Aj is centered on edges of i and k dimension but cell-centered along j
  // Ak is centered on edges of i and j dimension but cell-centered along k
  // Aj_right(k,j,i) = Aj(k+1/2,j,i-1/2)    Aj_left(k,j,i) = Aj(k-1/2,j,i-1/2)
  // Ak_right(k,j,i) = Ak(k,j+1/2,i-1/2)    Ak_left(k,j,i) = Ak(k,j-1/2,i-1/2)

  // get dimensions of interface field
  int fc_mx = bfield.length_dim0();
  int fc_my = bfield.length_dim1();
  int fc_mz = bfield.length_dim2();

  EnzoArray<enzo_float> Ak_left, Ak_right,Aj_right, Aj_left;

  if (dim == 0){
    Aj_right.initialize_subarray(Aj, 1, fc_mz+1, 0, fc_my, 0, fc_mx);
    Ak_right.initialize_subarray(Ak, 0, fc_mz, 1, fc_my+1, 0, fc_mx);
  } else if (dim == 1){
    Aj_right.initialize_subarray(Aj, 0, fc_mz, 0, fc_my, 1, fc_mx+1);
    Ak_right.initialize_subarray(Ak, 1, fc_mz+1, 0, fc_my, 0, fc_mx);
  } else {
    Aj_right.initialize_subarray(Aj, 0, fc_mz, 1, fc_my+1, 0, fc_mx);
    Ak_right.initialize_subarray(Ak, 0, fc_mz, 0, fc_my, 0, fc_mx+1);
  }

  Aj_left.initialize_subarray(Aj,  0,   fc_mz, 0, fc_my, 0, fc_mx);
  Ak_left.initialize_subarray(Ak,  0,   fc_mz, 0, fc_my, 0, fc_mx);

  for (int iz=0;iz<fc_mz;iz++){
    for (int iy=0; iy<fc_my; iy++){
      for (int ix=0; ix<fc_mx; ix++){
	// Bi(k,j,i-1/2) =
	//    ( Ak(    k, j+1/2, i-1/2) - Ak(    k, j-1/2, i-1/2) )/dj -
	//    ( Aj(k+1/2,     j, i-1/2) - Aj(k+1/2,     j, i-1/2) )/dk
	bfield(iz,iy,ix) = ((Ak_right(iz,iy,ix) - Ak_left(iz,iy,ix))/dj -
			    (Aj_right(iz,iy,ix) - Aj_left(iz,iy,ix))/dk);
	
      }
    }
  }
}

// Components of the magnetic field from the vector potential
// Bx(k, j, i-1/2) =
//    ( Az(    k, j+1/2, i-1/2) - Az(    k, j-1/2, i-1/2) )/dy -
//    ( Ay(k+1/2,     j, i-1/2) - Ay(k+1/2,     j, i-1/2) )/dz
// By(k, j-1/2, i) =
//    ( Ax(k+1/2, j-1/2,     i) - Ax(k+1/2, j-1/2,     i) )/dz -
//    ( Az(    k, j-1/2, i+1/2) - Az(    k, j-1/2, i-1/2) )/dx
// Bz(k-1/2, j, i) =
//    ( Ay(k-1/2,     j, i+1/2) - Ay(k-1/2,     j, i-1/2) )/dx -
//    ( Ax(k-1/2, j+1/2,     i) - Ax(k-1/2, j-1/2,     i) )/dy
  
void setup_bfield(Block * block, VectorInitializer &a, MeshPos &pos,
		  int mx, int my, int mz)
{
  EnzoBlock * enzo_block = enzo::block(block);

  EnzoArray<enzo_float> bfieldi_x, bfieldi_y, bfieldi_z;
  initialize_field_array(block, bfieldi_x, "bfieldi_x");
  initialize_field_array(block, bfieldi_y, "bfieldi_y");
  initialize_field_array(block, bfieldi_z, "bfieldi_z");

  // cell widths
  enzo_float dx = enzo_block->CellWidth[0];
  enzo_float dy = enzo_block->CellWidth[1];
  enzo_float dz = enzo_block->CellWidth[2];

  // allocate corner-centered arrays for the magnetic vector potentials
  // Ax, Ay, and Az are always cell-centered along the x, y, and z dimensions,
  // respectively
  EnzoArray<enzo_float> Ax, Ay, Az;
  Ax.initialize(mz+1,my+1,mx);
  Ay.initialize(mz+1,my,mx+1);
  Az.initialize(mz,my+1,mx+1);
  
  // Compute the Magnetic Vector potential at all points on the grid
  for (int iz=0; iz<mz+1; iz++){
    for (int iy=0; iy<my+1; iy++){
      for (int ix=0; ix<mx+1; ix++){

	enzo_float temp_ax, temp_ay, temp_az;

	if (ix != mx){
	  a(pos.x(iz,iy,ix),  pos.y_face(iz,iy,ix),  pos.z_face(iz,iy,ix),
	    Ax(iz,iy,ix), temp_ay, temp_az);
	}
	if (iy != my){
	  a(pos.x_face(iz,iy,ix),  pos.y(iz,iy,ix),  pos.z_face(iz,iy,ix),
	    temp_ay, Ay(iz,iy,ix), temp_az);
	}
	if (iz != mz){
	  a(pos.x_face(iz,iy,ix),  pos.y_face(iz,iy,ix),  pos.z(iz,iy,ix),
	    temp_ax, temp_ay, Az(iz,iy,ix));
	}
      }
    }
  }

  // Compute the Interface B-fields
  bfieldi_helper_(bfieldi_x, Ay, Az, 0, dy,dz);
  bfieldi_helper_(bfieldi_y, Az, Ax, 1, dz,dx);
  bfieldi_helper_(bfieldi_z, Ax, Ay, 2, dx,dy);
}


void setup_fluid(Block *block, ScalarInitializer &density_init,
		 ScalarInitializer &total_energy_init, 
		 VectorInitializer &momentum_init,
		 MeshPos &pos, int mx, int my, int mz)
{
  EnzoArray<enzo_float> density, total_energy;
  initialize_field_array(block, density, "density");
  initialize_field_array(block, total_energy, "total_energy");

  std::vector<EnzoArray<enzo_float>> momentum(3);
  initialize_field_array(block, momentum[0], "momentum_x");
  initialize_field_array(block, momentum[1], "momentum_y");
  initialize_field_array(block, momentum[2], "momentum_z");

  for (int iz=0;iz<mz;iz++){
    for (int iy=0; iy<my; iy++){
      for (int ix=0; ix<mx; ix++){
	enzo_float x,y,z;
	x = pos.x(iz,iy,ix); y = pos.y(iz,iy,ix); z = pos.z(iz,iy,ix);

	density(iz,iy,ix) = density_init(x,y,z);
	total_energy(iz,iy,ix) = total_energy_init(x,y,z);

	momentum_init(x,y,z,momentum[0](iz,iy,ix), momentum[1](iz,iy,ix),
		      momentum[2](iz,iy,ix));
      }
    }
  }
}


//----------------------------------------------------------------------

EnzoInitialLinearWave::EnzoInitialLinearWave(int cycle, double time,
					     double alpha, double beta,
					     double gamma,
					     int wave_type) throw()
  : Initial (cycle,time),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma),
    wave_type_(wave_type)
{}

//----------------------------------------------------------------------

void EnzoInitialLinearWave::enforce_block(Block * block,
					  const Hierarchy * hierarchy) throw()
{
  // Set up the test problem
  // This will only work for the VLCT method for an adiabatic fluid
  // does NOT support AMR
  
  // (PPML initial conditions are much more complicated)
  // Maybe don't pre-calculate alpha and beta

  Rotation rot(alpha_,beta_);

  EnzoBlock * enzo_block = enzo::block(block);
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];

  // Initialize object to compute positions
  MeshPos pos(enzo_block);

  // based on the wave type, setup the above scalar and vector initializers
  //RotatedScalarInitializer density(&rot;
  //				   new LinearScalarInitializer(backgroud,
  //							       eigenvalue,
  // 							       amplitude,
  //							       lambda));
  //RotatedVectorInitializer momentum(&rot;
  //				    new LinearVectorInitializer(back0, back1,
  //								back2, ev0, ev1,
  //								ev2, amplitude,
  //								lambda));

  // setup_bfield
  // setup_fluid
}
