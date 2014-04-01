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


// Pointers to dimension-independent functions
Number (*exact_solution_u)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_u)(const Point& p,
  const Parameters&, const std::string&, const std::string&);

Number (*exact_solution_x)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_x)(const Point& p,
  const Parameters&, const std::string&, const std::string&);

// Pointers to dimension-independent functions
Number (*exact_solution_p)(const Point& p,
  const Parameters&, const std::string&, const std::string&);
Gradient (*exact_derivative_p)(const Point& p,
  const Parameters&, const std::string&, const std::string&);


// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);
   
  std::string result_file_name ("data/cylinder_5567sym_100NT_10T_1_STAB");

  unsigned int n_timesteps = 100;
  Real time     = 0;
  Real end_time     = 10;

#if !THREED
  unsigned int N_eles=10;
#endif
#if THREED
  unsigned int N_eles=5;
#endif

if (argc >2){
read_options(n_timesteps,N_eles,result_file_name,argc, argv);
}

Real dt = end_time/n_timesteps;
#if THREED
  // Create a mesh.
  Mesh mesh(3);
#endif


/*
#if !THREED
    Mesh mesh;
  MeshTools::Generation::build_square (mesh,
                                       N_eles, N_eles,
                                       0., 1.,
                                       0., 1.,
                                       MESH_ELEMENT);
#endif
#if THREED
  MeshTools::Generation::build_cube (mesh,
                                       N_eles, N_eles, N_eles,
                                       0., 1.,
                                       0., 1.,
                                       0., 1.,
                                       HEX27);
#endif
*/


 // std::string mesh_file_name ("cylinder_sym728.msh");
  std::string mesh_file_name ("cylinder_5567sym.msh");
//   std::string mesh_file_name ("cylinder_sym737.msh");

  std::cout << mesh_file_name << std::endl;
  GmshIO(mesh).read(mesh_file_name);
  mesh.prepare_for_use();
  #if !PRES_STAB
  mesh.all_second_order(); //Need fpr P2
  #endif
  mesh.print_info();



  EquationSystems equation_systems (mesh);
  equation_systems.parameters.set<Real> ("dt")   = dt;

  ExodusII_IO exo= ExodusII_IO(mesh);

#if WRITE_TEC
TecplotIO tec= TecplotIO(equation_systems.get_mesh());
#endif
  
  TransientLinearImplicitSystem & system = 
    equation_systems.add_system<TransientLinearImplicitSystem> ("Last_non_linear_soln");

  system.add_variable ("s_u", DISP_ORDER,ELEMENT_TYPE);
  system.add_variable ("s_v", DISP_ORDER,ELEMENT_TYPE);
  #if THREED
  system.add_variable ("s_w", DISP_ORDER,ELEMENT_TYPE);
  #endif
  system.add_variable ("s_p", PRES_ORDER,ELEMENT_TYPE_PRESS);
  system.add_variable ("x", VEL_ORDER,ELEMENT_TYPE);
  system.add_variable ("y", VEL_ORDER,ELEMENT_TYPE);
  #if THREED
  system.add_variable ("z", VEL_ORDER,ELEMENT_TYPE);
  #endif

  system.attach_assemble_function (assemble_stokes);
  

  TransientLinearImplicitSystem & result =   equation_systems.add_system<TransientLinearImplicitSystem> ("result");
  result.add_variable ("u", DISP_ORDER,ELEMENT_TYPE);
  result.add_variable ("v", DISP_ORDER,ELEMENT_TYPE);
  #if THREED
  result.add_variable ("w", DISP_ORDER,ELEMENT_TYPE);
  #endif

  TransientLinearImplicitSystem & reference =   equation_systems.add_system<TransientLinearImplicitSystem> ("reference");
  reference.add_variable ("ref_u", DISP_ORDER,ELEMENT_TYPE);
  reference.add_variable ("ref_v", DISP_ORDER,ELEMENT_TYPE);
  #if THREED
  reference.add_variable ("ref_w", DISP_ORDER,ELEMENT_TYPE);
  #endif
  equation_systems.init ();
equation_systems.print_info();

  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 2500;
  equation_systems.parameters.set<Real>        ("linear solver tolerance") = TOLERANCE;
  
#if PETSC_MUMPS
  PetscLinearSolver<Number>* petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(system.get_linear_solver());
  PC pc = petsc_linear_solver->pc();
  int ierr = PCSetType(pc, PC_TYPE);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
  ierr = PCFactorSetMatSolverPackage(pc,SOLVER_NAME);
  CHKERRABORT(libMesh::COMM_WORLD,ierr);
#endif

//Set up exact solution
#if ANAL_2D
exact_solution_u = &exact_2D_solution_u;
exact_derivative_u = &exact_2D_derivative_u;

exact_solution_x = &exact_2D_solution_x;
exact_derivative_x = &exact_2D_derivative_x;

exact_solution_p = &exact_2D_solution_p;
exact_derivative_p = &exact_2D_derivative_p;

ExactSolution exact_sol_u(equation_systems);
exact_sol_u.attach_exact_value(exact_solution_u);
exact_sol_u.attach_exact_deriv(exact_derivative_u);

ExactSolution exact_sol_x(equation_systems);
exact_sol_x.attach_exact_value(exact_solution_x);
exact_sol_x.attach_exact_deriv(exact_derivative_x);

ExactSolution exact_sol_p(equation_systems);
exact_sol_p.attach_exact_value(exact_solution_p);
exact_sol_p.attach_exact_deriv(exact_derivative_p);

std::string u ("s_u");
std::string v ("s_v");
std::string p ("s_p");
std::string x ("x");
std::string y ("y");
#endif

  for (unsigned int t_step=1; t_step<=n_timesteps; ++t_step)
  {

    time += dt;
    equation_systems.parameters.set<Real> ("time") = time;
    double progress = (t_step+0.000000001) / (n_timesteps+0.000000001);
    equation_systems.parameters.set<Real>("progress") = progress;
    equation_systems.parameters.set<unsigned int>("step") = t_step; 

    std::cout << "\n\n*** Solving time step " << t_step << ", time = " << time <<  ", progress = " << progress << " ***" << std::endl;

    *system.old_local_solution = *system.current_local_solution;

  #if PETSC_MUMPS
      petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(system.get_linear_solver());
      pc = petsc_linear_solver->pc();
      ierr = PCSetType(pc, PC_TYPE);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
      ierr = PCFactorSetMatSolverPackage(pc,SOLVER_NAME);
      CHKERRABORT(libMesh::COMM_WORLD,ierr);
    #endif

    equation_systems.get_system("Last_non_linear_soln").solve();


    // How many iterations were required to solve the linear system?
    std::cout<<"Number of iterations: "<<system.n_linear_iterations()<<std::endl;        
    // What was the final residual of the linear system?
    std::cout<<"Residual: "<< system.final_linear_residual()<<std::endl;

    std::cout<<"About to update mesh (only on clpc59)" <<std::endl;

    //Update the mesh position - For cylinder only works on clpc59 !!!
    Mesh::node_iterator it_node = mesh.nodes_begin();
    const Mesh::node_iterator it_last_node = mesh.nodes_end();
    for ( ; it_node != it_last_node ; ++it_node)
    {
      Node* node_temp = *it_node;
      #if THREED
      for (unsigned int d = 0; d < 3; ++d) {
      #endif
      #if !THREED
      for (unsigned int d = 0; d < 2; ++d) {
      #endif
        unsigned int dest_dof = node_temp->dof_number(result.number(), d, 0);
        unsigned int dest_dof_system = node_temp->dof_number(result.number(), d, 0);
        Real value = (*node_temp)(d) + system.current_solution(dest_dof_system);
        result.current_local_solution->set(dest_dof, value);
        result.solution->set(dest_dof, value);

        reference.current_local_solution->set(dest_dof, (*node_temp)(d));
        reference.solution->set(dest_dof, (*node_temp)(d));
      }
    }
    
     #ifdef LIBMESH_HAVE_EXODUS_API
    //ExodusII_IO(mesh).write_equation_systems (result_file_name,                   
    std::stringstream file_name;
    file_name << result_file_name;
    file_name << std::setw(2) << std::setfill('0') << t_step;
    file_name << ".e-s.";
    file_name << std::setw(3) << std::setfill('0') << t_step+1;
    exo.write_timestep(file_name.str(), equation_systems,t_step+1,time); 
    std::cout<<"exodus "<< file_name.str() <<std::endl;

	exo.write_element_data(equation_systems);

    #endif 



 #if WRITE_TEC
  std::stringstream file_name_tec;
  file_name_tec << result_file_name << "_"<< t_step<< ".tec" ;
  tec.write_equation_systems (file_name_tec.str(),equation_systems);
std::cout<<"Wrote "<< file_name_tec.str() <<std::endl;
  #endif


 }


  return 0;
}


#include "assemble.h"
#include "assemble_stokes.cpp"
#include "exact_functions.cpp"
#include "read_options.cpp"
#include "test.cpp"
//#include "assemble_error.cpp"
