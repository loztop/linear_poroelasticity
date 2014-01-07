#include "assemble.h"

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <time.h>
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

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "dense_submatrix.h"
#include "dense_subvector.h"

// The definition of a geometric element
#include "elem.h"
#include <fe_interface.h>

#include "assemble.h"

void assemble_error (Real& H1_semi_error_disp, Real& Hdiv_semi_error_vel, Real& L2_error_press, Real& L2_error_displacement, Real& L2_error_velocity, EquationSystems& es,
                      const std::string& system_name)
{

  PerfLog perf_log("Assemble");

  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert (system_name == "Last_non_linear_soln");
  
  const Real dt    = es.parameters.get<Real>("dt");
  const Real progress    = es.parameters.get<Real>("progress");
  const Real time    = es.parameters.get<Real>("time");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  
  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem> ("Last_non_linear_soln");

 //TransientLinearImplicitSystem & system =    es.get_system<TransientLinearImplicitSystem> ("Stokes");
  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = system.variable_number ("s_u");
  const unsigned int v_var = system.variable_number ("s_v");
  #if THREED
  const unsigned int w_var = system.variable_number ("s_w");
  #endif
  const unsigned int p_var = system.variable_number ("s_p");
  const unsigned int x_var = system.variable_number ("x");
  const unsigned int y_var = system.variable_number ("y");
  #if THREED
  const unsigned int z_var = system.variable_number ("z");
  #endif
  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_disp_type = system.variable_type(u_var);
  FEType fe_vel_type = system.variable_type(x_var);

  // Get the Finite Element type for "p".
  FEType fe_pres_type = system.variable_type(p_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  AutoPtr<FEBase> fe_disp  (FEBase::build(dim, fe_disp_type));
  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
    
  // Build a Finite Element object of the specified type for
  // the pressure variables.
  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  
  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  
  //QGauss qrule (dim, fe_vel_type.default_quadrature_order());

 // QGauss qrule (dim, SIXTEENTH);
    QGauss qrule (dim, EIGHTH);

// QGauss qrule (dim, CONSTANT);

  // Tell the finite element objects to use our quadrature rule.
  fe_disp->attach_quadrature_rule (&qrule);
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  
  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe_vel->get_JxW();
  
  const std::vector<Point>& q_point = fe_vel->get_xyz();
  const std::vector<Point>& q_point_p = fe_pres->get_xyz();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi = fe_disp->get_dphi();
  const std::vector<std::vector<Real> >& phi = fe_disp->get_phi();
  const std::vector<std::vector<RealGradient> >& f_dphi = fe_vel->get_dphi();
  const std::vector<std::vector<Real> >& f_phi = fe_vel->get_phi();

  // The element shape functions for the pressure variable
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".


  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  #if THREED
  std::vector<unsigned int> dof_indices_w;
  #endif
  std::vector<unsigned int> dof_indices_p;
  std::vector<unsigned int> dof_indices_x;
  std::vector<unsigned int> dof_indices_y;
  #if THREED
  std::vector<unsigned int> dof_indices_z;
  #endif
  
  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users later
  // modify this program to include refinement, we will be safe and
  // will only consider the active elements; hence we use a variant of
  // the \p active_elem_iterator.

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  

//vectors needed to store bcs. (wha happened to condense ?)
std::vector<  int > rows;
std::vector<  int > pressure_rows;
std::vector< Real > rows_values;
std::vector< Real > pressure_rows_values;

std::vector< int> stab_dofs_rows;
std::vector< int> stab_dofs_cols;
std::vector<Real> stab_dofs_vals;

std::vector<Real> error_disp(4);
std::vector<Real> error_vel(4);
std::vector<Real> error_p(4);

//std::cout<< error_vals.size() <<std::endl;

  for ( ; el != end_el; ++el)
    {    

       test(6);

      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;
      
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      #if THREED
      dof_map.dof_indices (elem, dof_indices_w, w_var);
      #endif
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_x, x_var);
      dof_map.dof_indices (elem, dof_indices_y, y_var);
      #if THREED
      dof_map.dof_indices (elem, dof_indices_z, z_var);
      #endif

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();
      const unsigned int n_x_dofs = dof_indices_x.size(); 
      const unsigned int n_y_dofs = dof_indices_y.size();
      
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe_disp->reinit  (elem);
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);


      // Now we will build the element matrix.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

  
        //Displacement errors
        Real u_h = 0.;
        RealGradient grad_u_h;
        Real v_h = 0.;
        RealGradient grad_v_h;
        #if THREED
        Real w_h = 0.;
        RealGradient grad_w_h;
        #endif
        for (unsigned int i=0; i<n_u_dofs; i++)
          {
               u_h      += phi[i][qp]*system.current_solution  (dof_indices_u[i]);
               grad_u_h += dphi[i][qp]*system.current_solution (dof_indices_u[i]);
               v_h      += phi[i][qp]*system.current_solution  (dof_indices_v[i]);
               grad_v_h += dphi[i][qp]*system.current_solution (dof_indices_v[i]);   
              #if THREED
               w_h      += phi[i][qp]*system.current_solution  (dof_indices_w[i]);
               grad_w_h += dphi[i][qp]*system.current_solution (dof_indices_w[i]);
              #endif
          }

        Real exact_u = exact_2D_solution_u(q_point[qp], es.parameters,"null","void");
        Real error_u_sq =  (pow(u_h - exact_u,2));
        Real exact_v = exact_2D_solution_v(q_point[qp], es.parameters,"null","void");
        Real error_v_sq = (pow(v_h - exact_v,2));
        error_disp[0] += JxW[qp]*error_u_sq+JxW[qp]*error_v_sq;
        #if THREED
        Real exact_w = exact_2D_solution_w(q_point[qp], es.parameters,"null","void");
        Real error_w_sq = (pow(w_h - exact_w,2));
        error_disp[0] += JxW[qp]*error_w_sq;
        #endif

        RealGradient exact_gradu =exact_2D_derivative_u(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_u = grad_u_h - exact_gradu;
        RealGradient exact_gradv =exact_2D_derivative_v(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_v = grad_v_h - exact_gradv;
        error_disp[1] += JxW[qp]*error_grad_u.size_sq()+JxW[qp]*error_grad_v.size_sq();
        #if THREED
        RealGradient exact_gradw =exact_2D_derivative_w(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_w = grad_w_h - exact_gradw;
        error_disp[1] += JxW[qp]*error_grad_w.size_sq();
        #endif

      
        //Velocity errors
        Real x_h = 0.;
        RealGradient grad_x_h;
        Real y_h = 0.;
        RealGradient grad_y_h;
        #if THREED
        Real z_h = 0.;
        RealGradient grad_z_h;
        #endif
        for (unsigned int i=0; i<n_x_dofs; i++)
          {
               x_h      += phi[i][qp]*system.current_solution  (dof_indices_x[i]);
               grad_x_h += dphi[i][qp]*system.current_solution (dof_indices_x[i]);
               y_h      += phi[i][qp]*system.current_solution  (dof_indices_y[i]);
               grad_y_h += dphi[i][qp]*system.current_solution (dof_indices_y[i]);
               #if THREED
               z_h      += phi[i][qp]*system.current_solution  (dof_indices_z[i]);
               grad_z_h += dphi[i][qp]*system.current_solution (dof_indices_z[i]);
              #endif
          }

        Real exact_x = exact_2D_solution_x(q_point[qp], es.parameters,"null","void");
        Real error_x_sq = (pow(x_h - exact_x,2));
        Real exact_y = exact_2D_solution_y(q_point[qp], es.parameters,"null","void");
        Real error_y_sq = (pow(y_h - exact_y,2));
        error_vel[0] += JxW[qp]*error_x_sq+JxW[qp]*error_y_sq;
        #if THREED
        Real exact_z = exact_2D_solution_z(q_point[qp], es.parameters,"null","void");
        Real error_z_sq = (pow(z_h - exact_z,2));
        error_vel[0] += JxW[qp]*error_z_sq;
        #endif

        RealGradient exact_gradx =exact_2D_derivative_x(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_x = grad_x_h - exact_gradx;
        RealGradient exact_grady =exact_2D_derivative_y(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_y = grad_y_h - exact_grady;
        error_vel[1] += JxW[qp]*error_grad_x.size_sq()+JxW[qp]*error_grad_y.size_sq();
        #if THREED
        RealGradient exact_gradz =exact_2D_derivative_z(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_z = grad_z_h - exact_gradz;
        error_vel[1] += JxW[qp]*error_grad_z.size_sq();
        #endif

        RealGradient div_error_vel;
        #if !THREED
        div_error_vel(0) = grad_x_h(0) - exact_gradx(0)+grad_y_h(1) - exact_grady(1);
        #endif
        #if THREED
        div_error_vel(0) = grad_x_h(0) - exact_gradx(0)+grad_y_h(1) - exact_grady(1) +grad_z_h(2) - exact_gradz(2);
        #endif
        error_vel[2] += JxW[qp]*div_error_vel.size_sq();


        //Pressure errors
        Real p_h = 0.;
        for (unsigned int i=0; i<n_p_dofs; i++)
          {
               p_h      += psi[i][qp]*system.current_solution  (dof_indices_p[i]);
          }


        Real exact_p = exact_2D_solution_p(q_point_p[qp], es.parameters,"null","void");
        Real error_p_sq = (pow(p_h - exact_p,2));
        error_p[0] += JxW[qp]*error_p_sq;

} // end qp



} // end of element loop
  
  H1_semi_error_disp=sqrt(error_disp[1]);
  Hdiv_semi_error_vel=sqrt(error_vel[2]);  
  L2_error_press=sqrt(error_p[0]);
  L2_error_displacement=sqrt(error_disp[0]);
  L2_error_velocity=sqrt(error_vel[0]);

 // std::cout<<"H1_semi_error_disp "<< H1_semi_error_disp<<std::endl;
 // std::cout<<"Hdiv_semi_error_vel "<<Hdiv_semi_error_vel<<std::endl;
// std::cout<<"L2_error_disp "<<sqrt(error_disp[0])<<std::endl;
 //  std::cout<<"L2_error_press "<< sqrt(error_p[0])<<std::endl;

  // That's it.
  return;
}


