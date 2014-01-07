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

void assemble_error (EquationSystems& es,
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
  const unsigned int p_var = system.variable_number ("s_p");
  const unsigned int x_var = system.variable_number ("x");
  const unsigned int y_var = system.variable_number ("y");
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
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

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
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseMatrix<Number> Kstab;


  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke), Kux(Ke), Kuy(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke), Kvx(Ke), Kvy(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke), Kpx(Ke), Kpy(Ke),
    Kxu(Ke), Kxv(Ke), Kxp(Ke), Kxx(Ke), Kxy(Ke),
    Kyu(Ke), Kyv(Ke), Kyp(Ke), Kyx(Ke), Kyy(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe),
    Fx(Fe),
    Fy(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_p;
  std::vector<unsigned int> dof_indices_x;
  std::vector<unsigned int> dof_indices_y;
  
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
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_x, x_var);
      dof_map.dof_indices (elem, dof_indices_y, y_var);

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

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
      Kux.reposition (u_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs , n_u_dofs, n_x_dofs);
      Kuy.reposition (u_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_u_dofs, n_y_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
      Kvx.reposition (v_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs , n_v_dofs, n_x_dofs);
      Kvy.reposition (v_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_v_dofs, n_y_dofs);

      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
      Kpx.reposition (p_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs , n_p_dofs, n_x_dofs);
      Kpy.reposition (p_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_p_dofs, n_y_dofs);

      Kxu.reposition (p_var*n_u_dofs + n_p_dofs, u_var*n_u_dofs, n_x_dofs, n_u_dofs);
      Kxv.reposition (p_var*n_u_dofs + n_p_dofs, v_var*n_u_dofs, n_x_dofs, n_v_dofs);
      Kxp.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs, n_x_dofs, n_p_dofs);
      Kxx.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs , n_x_dofs, n_x_dofs);
      Kxy.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_x_dofs, n_y_dofs);

      Kyu.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, u_var*n_u_dofs, n_y_dofs, n_u_dofs);
      Kyv.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, v_var*n_u_dofs, n_y_dofs, n_v_dofs);
      Kyp.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs, n_y_dofs, n_p_dofs);
      Kyx.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs + n_p_dofs , n_y_dofs, n_x_dofs);
      Kyy.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_y_dofs, n_y_dofs);


      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);
      Fx.reposition (p_var*n_u_dofs + n_p_dofs, n_x_dofs);
      Fy.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, n_y_dofs);

    
  


      // Now we will build the element matrix.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
        test(7);
      

        //Displacement errors
        Real u_h = 0.;
        RealGradient grad_u_h;
        Real v_h = 0.;
        RealGradient grad_v_h;
        for (unsigned int i=0; i<n_u_dofs; i++)
          {
               u_h      += phi[i][qp]*system.current_solution  (dof_indices_u[i]);
               grad_u_h += dphi[i][qp]*system.current_solution (dof_indices_u[i]);
               v_h      += phi[i][qp]*system.current_solution  (dof_indices_v[i]);
               grad_v_h += dphi[i][qp]*system.current_solution (dof_indices_v[i]);
          }

        Real exact_u = exact_2D_solution_u(q_point[qp], es.parameters,"null","void");
        Real error_u_sq =  sqrt(pow(u_h - exact_u,2));
        
        Real exact_v = exact_2D_solution_v(q_point[qp], es.parameters,"null","void");
        
        Real error_v_sq = sqrt(pow(v_h - exact_v,2));
       // error_disp[0] += JxW[qp]*error_u_sq+JxW[qp]*error_v_sq;
        error_disp[0] += JxW[qp]*error_u_sq;

        RealGradient exact_gradu =exact_2D_derivative_u(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_u = grad_u_h - exact_gradu;
        RealGradient exact_gradv =exact_2D_derivative_v(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_v = grad_v_h - exact_gradv;
        error_disp[1] += JxW[qp]*error_grad_u.size_sq()+JxW[qp]*error_grad_v.size_sq();

        Real error_ux_sq = sqrt(pow(grad_u_h(0) - exact_gradu(0),2));
        Real error_vy_sq = sqrt(pow(grad_v_h(1) - exact_gradv(1),2));
        error_disp[2] += JxW[qp]*error_ux_sq+JxW[qp]*error_vy_sq;

        //Velocity errors
        Real x_h = 0.;
        RealGradient grad_x_h;
        Real y_h = 0.;
        RealGradient grad_y_h;
        for (unsigned int i=0; i<n_x_dofs; i++)
          {
               x_h      += phi[i][qp]*system.current_solution  (dof_indices_x[i]);
               grad_x_h += dphi[i][qp]*system.current_solution (dof_indices_x[i]);
               y_h      += phi[i][qp]*system.current_solution  (dof_indices_y[i]);
               grad_y_h += dphi[i][qp]*system.current_solution (dof_indices_y[i]);
          }

        Real exact_x = exact_2D_solution_x(q_point[qp], es.parameters,"null","void");
        Real error_x_sq = sqrt(pow(x_h - exact_x,2));
        Real exact_y = exact_2D_solution_y(q_point[qp], es.parameters,"null","void");
        Real error_y_sq = sqrt(pow(y_h - exact_y,2));
        //error_vel[0] += JxW[qp]*error_x_sq+JxW[qp]*error_y_sq;
        error_vel[0] += JxW[qp]*error_x_sq;

        RealGradient exact_gradx =exact_2D_derivative_x(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_x = grad_x_h - exact_gradx;
        RealGradient exact_grady =exact_2D_derivative_y(q_point[qp], es.parameters,"null","void");
        RealGradient error_grad_y = grad_y_h - exact_grady;
        //error_vel[1] += JxW[qp]*error_grad_x.size_sq()+JxW[qp]*error_grad_y.size_sq();
        error_vel[1] += JxW[qp]*error_grad_x.size_sq();

        Real error_xx_sq = sqrt(pow(grad_x_h(0) - exact_gradx(0),2));
        Real error_yy_sq = sqrt(pow(grad_y_h(1) - exact_grady(1),2));
        error_vel[2] += JxW[qp]*error_xx_sq+JxW[qp]*error_yy_sq;


        //Pressure errors
        Real p_h = 0.;
        for (unsigned int i=0; i<n_p_dofs; i++)
          {
               p_h      += psi[i][qp]*system.current_solution  (dof_indices_p[i]);
          }
        Real exact_p = exact_2D_solution_p(q_point_p[qp], es.parameters,"null","void");
        Real error_p_sq = sqrt(pow(p_h - exact_p,2));
        error_p[0] += JxW[qp]*error_p_sq;
} // end qp



} // end of element loop
  

    std::cout<<"disp L2 "<< sqrt(error_disp[0])<<std::endl;
    std::cout<<"disp H1 "<< sqrt(error_disp[0]+error_disp[1])<<std::endl;
    std::cout<<"disp div "<< sqrt(error_disp[0]+error_disp[2])<<std::endl;
    std::cout<<"vel L2 "<< sqrt(error_vel[0])<<std::endl;
    //std::cout<<"vel H1 "<< sqrt(error_vel[0]+error_vel[1])<<std::endl;
    std::cout<<"vel div "<< sqrt(error_vel[0]+error_vel[2])<<std::endl;
    std::cout<<"p L2 "<< sqrt(error_p[0])<<std::endl;

  // That's it.
  return;
}


