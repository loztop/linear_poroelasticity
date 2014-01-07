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
#include "dense_matrix_base.h"
#include "dense_vector_base.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "dense_submatrix.h"
#include "dense_subvector.h"
// The definition of a geometric element
#include "elem.h"
#include "assemble.h"

#define PT 1

void assemble_rhs (EquationSystems& es,
                      const std::string& system_name)
{

Real alpha=0.93;

#include "assemble_preamble.cpp"

  for ( ; el != end_el; ++el)
    {    

      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_x, x_var);
      dof_map.dof_indices (elem, dof_indices_y, y_var);
      #if THREED
      dof_map.dof_indices (elem, dof_indices_w, w_var);
      dof_map.dof_indices (elem, dof_indices_z, z_var);
      #endif

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();
      const unsigned int n_x_dofs = dof_indices_x.size(); 
      const unsigned int n_y_dofs = dof_indices_y.size();
      #if THREED
      const unsigned int n_w_dofs = dof_indices_w.size();
      const unsigned int n_z_dofs = dof_indices_z.size();
      #endif
      
      fe_disp->reinit  (elem);
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);


      Fe.resize (n_dofs);

     
      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);
      Fx.reposition (p_var*n_u_dofs + n_p_dofs, n_x_dofs);
      Fy.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, n_y_dofs);
      #if THREED
      Fw.reposition (w_var*n_u_dofs, n_w_dofs);
      Fz.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, n_y_dofs);
      #endif
    
      // Now we will build the element matrix.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {



        Point div_old_u;
        for (unsigned int l=0; l<n_u_dofs; l++)
        {
          div_old_u(0) += dphi[l][qp](0)*system.old_local_solution->el(dof_indices_u[l]);
          div_old_u(1) += dphi[l][qp](1)*system.old_local_solution->el(dof_indices_v[l]);
          #if THREED
          div_old_u(2) += dphi[l][qp](2)*system.old_local_solution->el(dof_indices_w[l]);
          #endif
        }


       /// std::cout<<div_old_u(0)<<std::endl;
		
		/*
        Real p_old=0;
        for (unsigned int l=0; l<n_p_dofs; l++)
        {
          p_old += psi[l][qp]*system.current_local_solution->el(dof_indices_p[l]);
        }
        Point z_old;
        for (unsigned int l=0; l<n_x_dofs; l++)
        {
          z_old(0) += f_phi[l][qp]*system.old_local_solution->el(dof_indices_x[l]);
          z_old(1) += f_phi[l][qp]*system.old_local_solution->el(dof_indices_y[l]);
        }
        Point grad_u_old;
        for (unsigned int l=0; l<n_u_dofs; l++)
        {
          grad_u_old(0) += dphi[l][qp](0)*system.old_local_solution->el(dof_indices_u[l]);
          grad_u_old(1) += dphi[l][qp](1)*system.old_local_solution->el(dof_indices_u[l]);
        }
        Point grad_v_old;
        for (unsigned int l=0; l<n_u_dofs; l++)
        {
          grad_v_old(0) += dphi[l][qp](0)*system.old_local_solution->el(dof_indices_v[l]);
          grad_v_old(1) += dphi[l][qp](1)*system.old_local_solution->el(dof_indices_v[l]);
        }
        */
        

		//add some gravity 
		for (unsigned int i=0; i<n_u_dofs; i++){
         // Fv(i) += -2*JxW[qp]*phi[i][qp];
        }

 
//Source term
          for (unsigned int i=0; i<n_p_dofs; i++){
            #if THREED
            Fp(i) += alpha*(div_old_u(0)+div_old_u(1)+div_old_u(2))*JxW[qp]*psi[i][qp];
            #endif
            #if !THREED
            Fp(i) +=  alpha*(div_old_u(0)+div_old_u(1))*JxW[qp]*psi[i][qp];
			
            #endif
            #if ANAL_2D
            Fp(i) +=  forcing_function_2D(q_point[qp], es.parameters)*JxW[qp]*psi[i][qp];
            #endif
          }
  
} // end qp


//Pressure jump stabilisation.
#if USE_STAB  
 	std::vector<unsigned int> stab_dofs_cols2;
	std::vector<Real> stab_dofs_vals2;


	std::vector<Real> stab_dofs_vals_rhs;

    for (unsigned int s=0; s<elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {   
				//Only do something on the interior edges.
      } //if (elem->neighbor(s) == NULL)
			else{

      const Elem* neighbor = elem->neighbor(s);

			std::vector<unsigned int> neighbor_dof_indices_p;
      dof_map.dof_indices(neighbor, neighbor_dof_indices_p,p_var);
      const unsigned int n_neighbor_dofs_p = neighbor_dof_indices_p.size();

			Real hmax=(*elem).hmax();
      Real hmin=(*elem).hmin();

			Real delta=pow(dt,delta_power)*PEN_STAB;
      #if !THREED
			Real factor=-delta*(hmax)*(hmax);;
      #endif
      
      #if THREED
      Real factor=-delta*(hmax*hmax*hmax);
      #endif

			//Accumulate RHS
			std::vector< Real > p_old;
			system.old_local_solution->get(dof_indices_p,p_old);
			std::vector< Real > p_old_n;
			system.old_local_solution->get(neighbor_dof_indices_p,p_old_n);
      
			stab_dofs_vals_rhs.push_back(-p_old[0]*factor);
			stab_dofs_vals_rhs.push_back(p_old_n[0]*factor);
      //perf_log.pop("push back");

			}
		}


//add rhs
DenseVector<Number> rhs_stab;
rhs_stab.resize(dof_indices_p.size());
 	
for (int i=0; i < stab_dofs_vals_rhs.size(); i++) {
	 	rhs_stab(0)+=stab_dofs_vals_rhs[i];
	}

std::vector<unsigned int> stab_dofs_rows_rhs;
stab_dofs_rows_rhs.push_back(dof_indices_p[0]);
system.rhs->add_vector(rhs_stab, stab_dofs_rows_rhs);
  #endif 

#if THREED
//  #include "assemble_stokes_bcs_p1p1p0_anal_sine.cpp"
//  #include "pin_pressure.cpp"
	
#endif

  #if !THREED

//  #include "assemble_stokes_bcs_p1p1p0_anal_sine.cpp"
//  #include "square_pin_pressure.cpp"
//  	#include "cube_squash_bcs.cpp"
#include "canteliver_bcs.cpp"
#include "canteliver_traction_bcs.cpp"
#endif

  system.rhs->add_vector    (Fe, dof_indices);


} // end of element loop
  
	//Apply BCS
    system.rhs->close();
    system.matrix->zero_rows(rows, 1.0);
    for (int i=0; i < rows.size(); i++) {
      system.rhs->set(rows[i],rows_values[i]);
    }

		system.matrix->close();
    system.rhs->close();

    std::cout<<pressure_rows.size()<<std::endl;
    system.matrix->zero_rows(pressure_rows, 1.0);
    for (int i=0; i < pressure_rows.size(); i++) {
      system.rhs->set(pressure_rows[i],pressure_rows_values[i]);
    }
    system.matrix->close();
    system.rhs->close();



    std::cout<<"Assemble rhs->l2_norm () "<<system.rhs->l2_norm ()<<std::endl;

  return;
}


