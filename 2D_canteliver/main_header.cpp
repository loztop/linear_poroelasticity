//iHat*u + jHat*v + kHat*w

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
#include <tecplot_io.h>
#include "error_vector.h"
#include "exact_solution.h"
#include <petsc_linear_solver.h>
#include "elem.h"
using namespace libMesh;
#include "gmsh_io.h"
#include "assemble.h"
#include <iostream>
#include <fstream>
using namespace std;



// Pointers to dimension-independent functions
Number (*exact_solution_u)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_u)(const Point& p,
  const Parameters&, const std::string&, const std::string&);

Number (*exact_solution_v)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_v)(const Point& p,
  const Parameters&, const std::string&, const std::string&);

#if THREED
Number (*exact_solution_w)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_w)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
#endif

Number (*exact_solution_x)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_x)(const Point& p,
  const Parameters&, const std::string&, const std::string&);

Number (*exact_solution_y)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_y)(const Point& p,
  const Parameters&, const std::string&, const std::string&);

#if THREED
Number (*exact_solution_z)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_z)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
#endif


// Pointers to dimension-independent functions
Number (*exact_solution_p)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_p)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
