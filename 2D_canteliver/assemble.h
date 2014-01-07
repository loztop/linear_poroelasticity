#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_

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



// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "dense_submatrix.h"
#include "dense_subvector.h"

// The definition of a geometric element
#include "elem.h"

#include "assemble.h"


#define THREED 0

#define SOLVER_NAME "mumps"
#define PC_TYPE PCLU
#define PETSC_MUMPS 1

#define WRITE_TEC 0
#define EXODUS 1


#define TIME 1


#define ANAL_2D 0

#define FULL_ELASTIC 1

#define NITSCHE 1

#define PEN_STAB 0.00005

//#define PEN_STAB 0.00000001

#define PEN_BC 100

#define E 1.0e+5
#define NU 0.4
#define KPERM 1.0e-7




#define DISP_ORDER FIRST
#define VEL_ORDER FIRST
#define PRES_ORDER CONSTANT
#define ELEMENT_TYPE LAGRANGE
#define ELEMENT_TYPE_PRESS MONOMIAL
#define MESH_ELEMENT TRI3
#define USE_STAB 1


/*
#define ANAL_2D 0
#define PRES_STAB 0
#define DISP_ORDER SECOND
#define VEL_ORDER SECOND
#define PRES_ORDER FIRST
#define ELEMENT_TYPE LAGRANGE
#define ELEMENT_TYPE_PRESS LAGRANGE
#define MESH_ELEMENT TRI6
*/

using namespace libMesh;

void assemble_stiffness (EquationSystems& es,
                      const std::string& system_name);

void assemble_rhs (EquationSystems& es,
                      const std::string& system_name);

void assemble_error(Real& H1_semi_error_disp, Real& Hdiv_semi_error_vel, Real& L2_error_press,Real& L2_error_displacement,Real& L2_error_velocity,  EquationSystems& es,
                      const std::string& system_name);

void read_options(unsigned int &  n_timesteps, unsigned int &  N_eles, std::string& output_file_name, std::string& result_file_name, int& argc, char**& argv) ;


void read_parameters(EquationSystems& es, int& argc, char**& argv) ;



void test(int a);

Number exact_2D_solution_u(const Point& p,  
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&,  // sys_name, not needed
                         const std::string&); // unk_name, not needed);

Number exact_2D_solution_v(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&,  // sys_name, not needed
                         const std::string&); // unk_name, not needed);

Number exact_2D_solution_w(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&,  // sys_name, not needed
                         const std::string&); // unk_name, not needed);

Gradient exact_2D_derivative_u(const Point& p,
  const Parameters& parameters, const std::string&, const std::string&);


Gradient exact_2D_derivative_v(const Point& p,
  const Parameters& parameters, const std::string&, const std::string&);

Gradient exact_2D_derivative_w(const Point& p,
  const Parameters& parameters, const std::string&, const std::string&);



Number exact_2D_solution_x(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&,  // sys_name, not needed
                         const std::string&); // unk_name, not needed);

Number exact_2D_solution_y(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&,  // sys_name, not needed
                         const std::string&); // unk_name, not needed);

Number exact_2D_solution_z(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&,  // sys_name, not needed
                         const std::string&); // unk_name, not needed);

Gradient exact_2D_derivative_x(const Point& p,
  const Parameters& parameters, const std::string&, const std::string&);


Gradient exact_2D_derivative_y(const Point& p,
  const Parameters& parameters, const std::string&, const std::string&);

Gradient exact_2D_derivative_z(const Point& p,
  const Parameters& parameters, const std::string&, const std::string&);



Number exact_2D_solution_p(const Point& p,
                         const Parameters& parameters,   // parameters, not needed
                         const std::string&,  // sys_name, not needed
                         const std::string&); // unk_name, not needed);

// Prototypes for calculation of the gradient of the exact solution.  
// Necessary for setting boundary conditions in H^2_0 and testing
// H^1 convergence of the solution
Gradient exact_2D_derivative_p(const Point& p,
  const Parameters& parameters, const std::string&, const std::string&);


Number forcing_function_2D(const Point& p, const Parameters& parameters);

#endif 
