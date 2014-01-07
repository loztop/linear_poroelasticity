#include "assemble.h"

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "exodusII_io.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "linear_implicit_system.h"
#include "transient_system.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"

// The definition of a geometric element
#include "elem.h"

#include "assemble.h"

#define FAC 1

#define POLY 0
#define SIN 1
#define PIE 3.1415926535897932384626433832

#define FIX_T 1

Number exact_2D_solution_u(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{ 

   Real t =  parameters.get<Real>("time");



  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  // analytic solution value
#if POLY
  return (20*x*y*y*y)/FAC;
#endif
#if SIN && ! THREED && !TIME
  return ((-1.0/(4.0*PIE))*cos(2.0*PIE*x)*sin(2.0*PIE*y))/FAC;
#endif


#if SIN && ! THREED && TIME
  return ((-1.0/(4.0*PIE))*cos(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*t))/FAC;
#endif

#if THREED && !TIME
  return ((-1.0/(6.0*PIE))*cos(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z))/FAC;
#endif


#if THREED && TIME
  return ((-1.0/(6.0*PIE))*cos(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
#endif

}

Number exact_2D_solution_v(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

  // analytic solution value


#if POLY
  return (5*x*x*x*x - 5*y*y*y*y)/FAC;
#endif
#if SIN && ! THREED && !TIME
  return ((-1.0/(4.0*PIE))*sin(2.0*PIE*x)*cos(2.0*PIE*y))/FAC;
#endif

#if SIN && ! THREED && TIME
  return ((-1.0/(4.0*PIE))*sin(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*t))/FAC;
#endif

#if THREED && !TIME
  return ((-1.0/(6.0*PIE))*sin(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*z))/FAC;
#endif


#if THREED && TIME
  return ((-1.0/(6.0*PIE))*sin(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
#endif

}

Number exact_2D_solution_w(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");


#if THREED && !TIME
  return ((-1.0/(6.0*PIE))*sin(2.0*PIE*x)*sin(2.0*PIE*y)*cos(2.0*PIE*z))/FAC;
#endif

#if THREED && TIME
  return ((-1.0/(6.0*PIE))*sin(2.0*PIE*x)*sin(2.0*PIE*y)*cos(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
#endif


}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_u(const Point& p,
                             const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");

  // First derivatives to be returned.
  Gradient gradu;

#if POLY
  gradu(0) = (20*y*y*y)/FAC;
  gradu(1) = (60*x*y*y)/FAC;
  return gradu;
#endif

#if SIN && !THREED && !TIME
  gradu(0) = ((1.0/2.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y))/FAC;
  gradu(1) = ((-1.0/2.0)*cos(2.0*PIE*x)*cos(2.0*PIE*y))/FAC;
  return gradu;
#endif

#if SIN && ! THREED && TIME
  gradu(0) = ((1.0/2.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*t))/FAC;
  gradu(1) = ((-1.0/2.0)*cos(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*t))/FAC;
  return gradu;
#endif

#if THREED && !TIME
  gradu(0) = ((1.0/3.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z))/FAC;
  gradu(1) = ((-1.0/3.0)*cos(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*z))/FAC;
  gradu(2) = ((-1.0/3.0)*cos(2.0*PIE*x)*sin(2.0*PIE*y)*cos(2.0*PIE*z))/FAC;
  return gradu;
#endif

#if THREED && TIME
  gradu(0) = ((1.0/3.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  gradu(1) = ((-1.0/3.0)*cos(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  gradu(2) = ((-1.0/3.0)*cos(2.0*PIE*x)*sin(2.0*PIE*y)*cos(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  return gradu;
#endif

}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_v(const Point& p,
                             const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");

  // First derivatives to be returned.
  Gradient gradv;

#if POLY
  gradv(0) = (20*x*x*x)/FAC;
  gradv(1) = (-20*y*y*y)/FAC;

  return gradv;
#endif

#if SIN && !THREED && !TIME
  gradv(0) = ((-1.0/2.0)*cos(2.0*PIE*x)*cos(2.0*PIE*y))/FAC;
  gradv(1) = ((1.0/2.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y))/FAC;
  return gradv;
#endif

#if SIN && !THREED && TIME
  gradv(0) = ((-1.0/2.0)*cos(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*t))/FAC;
  gradv(1) = ((1.0/2.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*t))/FAC;
  return gradv;
#endif

#if THREED && !TIME
  gradv(0) = ((-1.0/3.0)*cos(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*z))/FAC;
  gradv(1) = ((1.0/3.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z))/FAC;
  gradv(2) = ((-1.0/3.0)*sin(2.0*PIE*x)*cos(2.0*PIE*y)*cos(2.0*PIE*z))/FAC;
  return gradv;
#endif

#if THREED && TIME
  gradv(0) = ((-1.0/3.0)*cos(2.0*PIE*x)*cos(2.0*PIE*y)*sin(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  gradv(1) = ((1.0/3.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  gradv(2) = ((-1.0/3.0)*sin(2.0*PIE*x)*cos(2.0*PIE*y)*cos(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  return gradv;
#endif
}

// We now define the gradient of the exact solution
Gradient exact_2D_derivative_w(const Point& p,
                             const Parameters& parameters,   // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

  // First derivatives to be returned.
  Gradient gradw;

#if THREED && !TIME
  gradw(0) = ((-1.0/3.0)*cos(2.0*PIE*x)*sin(2.0*PIE*y)*cos(2.0*PIE*z))/FAC;
  gradw(1) = ((-1.0/3.0)*sin(2.0*PIE*x)*cos(2.0*PIE*y)*cos(2.0*PIE*z))/FAC;
  gradw(2) = ((1.0/3.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z))/FAC;
  return gradw;
#endif

#if THREED && TIME
  gradw(0) = ((-1.0/3.0)*cos(2.0*PIE*x)*sin(2.0*PIE*y)*cos(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  gradw(1) = ((-1.0/3.0)*sin(2.0*PIE*x)*cos(2.0*PIE*y)*cos(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  gradw(2) = ((1.0/3.0)*sin(2.0*PIE*x)*sin(2.0*PIE*y)*sin(2.0*PIE*z)*sin(2.0*PIE*t))/FAC;
  return gradw;
#endif
}



Number exact_2D_solution_x(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");


  // analytic solution value
#if POLY
  return (-120*x*y)/FAC;
#endif
#if SIN && !THREED && !TIME
  return (-KPERM*2*PIE*cos(2*PIE*x)*sin(2*PIE*y))/FAC;
#endif

#if SIN && !THREED && TIME
  return (-KPERM*2*PIE*cos(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*t))/FAC;
#endif

#if THREED  && !TIME
  return (-KPERM*2*PIE*cos(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC;
#endif

#if THREED  && TIME
  return (-KPERM*2*PIE*cos(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
#endif
}

Number exact_2D_solution_y(const Point& p ,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{

  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");


  // analytic solution value

#if POLY
  return (60*y*y -60*x*x)/FAC;
#endif

#if SIN && !THREED && !TIME
  return (-KPERM*2*PIE*sin(2*PIE*x)*cos(2*PIE*y))/FAC;
#endif

#if (SIN && !THREED) && TIME
  return (-KPERM*2*PIE*sin(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*t))/FAC;

  
#endif

#if THREED && !TIME
  return (-KPERM*2*PIE*sin(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*z))/FAC;
#endif

#if THREED && TIME
  return (-KPERM*2*PIE*sin(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
#endif
}

Number exact_2D_solution_z(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

  // analytic solution value
#if THREED && !TIME
  return (-KPERM*2*PIE*sin(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z))/FAC;
#endif

#if THREED && TIME
  return (-KPERM*2*PIE*sin(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z)*sin(2*PIE*t))/FAC;
#endif
}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_x(const Point& p,
                             const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
   Real t =  parameters.get<Real>("time");

  // First derivatives to be returned.
  Gradient gradx;

#if POLY
  gradx(0) = (120*y)/FAC;
  gradx(1) = (120*x)/FAC;

  return gradx;
#endif

#if SIN && !THREED && !TIME
  gradx(0) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y))/FAC;
  gradx(1) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*cos(2*PIE*y))/FAC;
  return gradx;
#endif

#if SIN && !THREED && TIME
  gradx(0) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*t))/FAC;
  gradx(1) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*t))/FAC;
  return gradx;
#endif

#if THREED && !TIME
  gradx(0) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC;
  gradx(1) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*z))/FAC;
  gradx(2) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z))/FAC;
  return gradx;
#endif

#if THREED && TIME
  gradx(0) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
  gradx(1) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
  gradx(2) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z)*sin(2*PIE*t))/FAC;
  return gradx;
#endif
}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_y(const Point& p,
                             const Parameters& parameters,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");


  // First derivatives to be returned.
  Gradient grady;

#if POLY
  grady(0) = (120*x)/FAC;
  grady(1) = (-120*y)/FAC;

  return grady;
#endif

#if SIN && !THREED && !TIME
  grady(0) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*cos(2*PIE*y))/FAC;
  grady(1) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y))/FAC;
  return grady;
#endif

#if SIN && !THREED && TIME
  grady(0) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*t))/FAC;
  grady(1) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*t))/FAC;
  return grady;
#endif

#if THREED && !TIME
  grady(0) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*z))/FAC;
  grady(1) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC;
  grady(2) = (-KPERM*4*PIE*PIE*sin(2*PIE*x)*cos(2*PIE*y)*cos(2*PIE*z))/FAC;
  return grady;
#endif


#if THREED && TIME
  grady(0) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
  grady(1) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
  grady(2) = (-KPERM*4*PIE*PIE*sin(2*PIE*x)*cos(2*PIE*y)*cos(2*PIE*z)*sin(2*PIE*t))/FAC;
  return grady;
#endif
}



// We now define the gradient of the exact solution
Gradient exact_2D_derivative_z(const Point& p,
                             const Parameters& parameters,   // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");

  // First derivatives to be returned.
  Gradient gradz;


#if THREED && !TIME
  gradz(0) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z))/FAC;
  gradz(1) = (-KPERM*4*PIE*PIE*sin(2*PIE*x)*cos(2*PIE*y)*cos(2*PIE*z))/FAC;
  gradz(2) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC;
  return gradz;
#endif


#if THREED && TIME
  gradz(0) = (-KPERM*4*PIE*PIE*cos(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z)*sin(2*PIE*t))/FAC;
  gradz(1) = (-KPERM*4*PIE*PIE*sin(2*PIE*x)*cos(2*PIE*y)*cos(2*PIE*z)*sin(2*PIE*t))/FAC;
  gradz(2) = (KPERM*4*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
  return gradz;
#endif
}

Number exact_2D_solution_p(const Point& p,
                         const Parameters& parameters,  // parameters, not needed
                         const std::string&, // sys_name, not needed
                         const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);
  Real t =  parameters.get<Real>("time");


#if POLY
  return (60.*x*x*y - 20*y*y*y)/FAC;
#endif
#if SIN && !THREED && !TIME
  return (sin(2*PIE*x)*sin(2*PIE*y))/FAC;
#endif

#if SIN && !THREED && TIME
  
  return (sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*t))/FAC;

#endif

#if THREED && !TIME
  return (sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC;
#endif

#if THREED && TIME
  return (sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC;
#endif

}


// We now define the gradient of the exact solution
Gradient exact_2D_derivative_p(const Point& p,
                             const Parameters&,  // parameters, not needed
                             const std::string&, // sys_name, not needed
                             const std::string&) // unk_name, not needed
{
  // x and y coordinates in space
  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  // First derivatives to be returned.
  Gradient gradp;


#if POLY
  gradp(0) = (120*x*y)/FAC;
  gradp(1) = (60*x*x-60*y*y)/FAC;

  return gradp;
#endif

#if SIN && !THREED
  gradp(0) = (-2*PIE*cos(2*PIE*x)*sin(2*PIE*y))/FAC;
  gradp(1) = (-2*PIE*sin(2*PIE*x)*cos(2*PIE*y))/FAC;
  return gradp;
#endif

#if THREED
  gradp(0) = (-2*PIE*cos(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC;
  gradp(1) = (-2*PIE*sin(2*PIE*x)*cos(2*PIE*y)*sin(2*PIE*z))/FAC;
  gradp(2) = (-2*PIE*sin(2*PIE*x)*sin(2*PIE*y)*cos(2*PIE*z))/FAC;
  return gradp;
#endif

}

Number forcing_function_2D(const Point& p, const Parameters& parameters)
{
    const Real x = p(0);
    const Real y = p(1);
    const Real z = p(2);
   Real t =  parameters.get<Real>("time");
  const Real dt =  parameters.get<Real>("dt");


  #if POLY
    return 0;
  #endif
  #if SIN  && ! THREED && !TIME

//BOTH STOKES AND DARCY DIV terms
  return (sin(2*PIE*x)*sin(2*PIE*y))/FAC   + (8*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y))/FAC   ;
  
#endif


 #if SIN  && ! THREED && TIME

  return (dt*sin(2*PIE*x)*sin(2*PIE*y)*2*PIE*cos(2*PIE*t))/FAC   + (dt*KPERM*8*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*t))/FAC   ;
   
     #endif


   #if THREED && !TIME
  //BOTH STOKES AND DARCY DIV terms
  return (sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC + (KPERM*12*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z))/FAC   ;
  #endif


  #if THREED && TIME
  //BOTH STOKES AND DARCY DIV terms
  return (dt*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z)*2*PIE*cos(2*PIE*t))/FAC + (dt*KPERM*12*PIE*PIE*sin(2*PIE*x)*sin(2*PIE*y)*sin(2*PIE*z)*sin(2*PIE*t))/FAC   ;

  #endif

  }
