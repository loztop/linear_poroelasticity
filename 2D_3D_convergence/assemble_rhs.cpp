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
#include "assemble_preamble.cpp"

/*

  
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



 
//Source term
          for (unsigned int i=0; i<n_p_dofs; i++){
            #if THREED
            Fp(i) += (div_old_u(0)+div_old_u(1)+div_old_u(2))*JxW[qp]*psi[i][qp];
            #endif
            #if !THREED
            Fp(i) +=  (div_old_u(0)+div_old_u(1))*JxW[qp]*psi[i][qp];
			
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
	  Real delta=dt*DELTA;

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
  #include "assemble_stokes_bcs_p1p1p0_anal_sine.cpp"
  #include "pin_pressure.cpp"
	
#endif

  #if !THREED

  #include "assemble_stokes_bcs_p1p1p0_anal_sine.cpp"
  #include "square_pin_pressure.cpp"

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
*/



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
      #if THREED
      Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
      Kuz.reposition (u_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_u_dofs, n_z_dofs);
      #endif

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
      Kvx.reposition (v_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs , n_v_dofs, n_x_dofs);
      Kvy.reposition (v_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_v_dofs, n_y_dofs);
      #if THREED
      Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);
      Kuz.reposition (v_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_u_dofs, n_z_dofs);
      #endif

      #if THREED
      Kwu.reposition (w_var*n_w_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kwv.reposition (w_var*n_w_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kwp.reposition (w_var*n_w_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
      Kwx.reposition (w_var*n_w_dofs, p_var*n_u_dofs + n_p_dofs , n_v_dofs, n_x_dofs);
      Kwy.reposition (w_var*n_w_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_v_dofs, n_y_dofs);
      Kww.reposition (w_var*n_w_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);
      Kwz.reposition (w_var*n_w_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_u_dofs, n_z_dofs);
      #endif

      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
      Kpx.reposition (p_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs , n_p_dofs, n_x_dofs);
      Kpy.reposition (p_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_p_dofs, n_y_dofs);
      #if THREED
      Kpw.reposition (p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_w_dofs);
      Kpz.reposition (p_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_p_dofs, n_z_dofs);
      #endif

      Kxu.reposition (p_var*n_u_dofs + n_p_dofs, u_var*n_u_dofs, n_x_dofs, n_u_dofs);
      Kxv.reposition (p_var*n_u_dofs + n_p_dofs, v_var*n_u_dofs, n_x_dofs, n_v_dofs);
      Kxp.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs, n_x_dofs, n_p_dofs);
      Kxx.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs , n_x_dofs, n_x_dofs);
      Kxy.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_x_dofs, n_y_dofs);
      #if THREED
      Kxw.reposition (p_var*n_u_dofs + n_p_dofs, w_var*n_u_dofs, n_x_dofs, n_w_dofs);
      Kxz.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_x_dofs, n_z_dofs);
      #endif


      Kyu.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, u_var*n_u_dofs, n_y_dofs, n_u_dofs);
      Kyv.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, v_var*n_u_dofs, n_y_dofs, n_v_dofs);
      Kyp.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs, n_y_dofs, n_p_dofs);
      Kyx.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs + n_p_dofs , n_y_dofs, n_x_dofs);
      Kyy.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_y_dofs, n_y_dofs);
      #if THREED
      Kyw.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, w_var*n_u_dofs, n_x_dofs, n_w_dofs);
      Kyz.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_x_dofs, n_z_dofs);
      #endif

      #if THREED
      Kzu.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, u_var*n_u_dofs, n_y_dofs, n_u_dofs);
      Kzv.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, v_var*n_u_dofs, n_y_dofs, n_v_dofs);
      Kzp.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, p_var*n_u_dofs, n_y_dofs, n_p_dofs);
      Kzx.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, p_var*n_u_dofs + n_p_dofs , n_y_dofs, n_x_dofs);
      Kzy.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_y_dofs, n_y_dofs);
      Kzw.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, w_var*n_u_dofs, n_x_dofs, n_w_dofs);
      Kzz.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_x_dofs, n_z_dofs);
      #endif



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
        test(7);


        Point div_old_u;
        for (unsigned int l=0; l<n_u_dofs; l++)
        {
          div_old_u(0) += dphi[l][qp](0)*system.old_local_solution->el(dof_indices_u[l]);
          div_old_u(1) += dphi[l][qp](1)*system.old_local_solution->el(dof_indices_v[l]);
          #if THREED
          div_old_u(2) += dphi[l][qp](2)*system.old_local_solution->el(dof_indices_w[l]);
          #endif
        }

            //Elastcity (Mixture) momentum equation




//Source term

          for (unsigned int i=0; i<n_p_dofs; i++){
            #if THREED
            Fp(i) += (div_old_u(0)+div_old_u(1)+div_old_u(2))*JxW[qp]*psi[i][qp];
            #endif

            #if !THREED
            Fp(i) += (div_old_u(0)+div_old_u(1))*JxW[qp]*psi[i][qp];
            #endif
          }

        #if ANAL_2D
					for (unsigned int i=0; i<n_p_dofs; i++){
            Fp(i) += forcing_function_2D(q_point[qp], es.parameters)*JxW[qp]*psi[i][qp];
          }
        #endif



  
} // end qp

	//Pressure jump stabilisation.
#if USE_STAB  
 	std::vector<unsigned int> stab_dofs_cols2;
	std::vector<Real> stab_dofs_vals2;
    for (unsigned int s=0; s<elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {   
				//Only do something on the interior edges.
      } //if (elem->neighbor(s) == NULL)
			else{

      const Elem* neighbor = elem->neighbor(s);
   		AutoPtr<Elem> side (elem->build_side(s));


      std::vector<Real> node_co_ords_x;
      std::vector<Real> node_co_ords_y;
      for (unsigned int n=0; n<elem->n_nodes(); n++)
      {
          Node *e_node = elem->get_node(n);
          for (unsigned int m=0; m<neighbor->n_nodes(); m++)
          {
            Node *n_node = neighbor->get_node(m);
            if( ( (*n_node)(0)==(*e_node)(0) ) && ( (*n_node)(1)==(*e_node)(1) ) ){
              node_co_ords_x.push_back((*n_node)(0));
              node_co_ords_y.push_back((*n_node)(1));
              }
          }
      }
     Real side_length=pow ( pow((node_co_ords_x[0]-node_co_ords_x[1]),2) + pow((node_co_ords_y[0]-node_co_ords_y[1]),2),0.5) ;

			std::vector<unsigned int> neighbor_dof_indices_p;
      dof_map.dof_indices(neighbor, neighbor_dof_indices_p,p_var);
      const unsigned int n_neighbor_dofs_p = neighbor_dof_indices_p.size();
			Real hmax=(*elem).hmax();
      Real hmin=(*elem).hmin();
			//Real vol=(*elem).volume();
  const Real DELTA    = es.parameters.get<Real>("DELTA");
  Real delta=dt*DELTA;

		//	Real delta=dt*PEN_STAB;
      #if !THREED
			Real factor=-delta*(hmax)*(hmax);;
      //Real factor=-delta*(side_length);
      #endif
      
      #if THREED
      Real factor=-delta*(hmax*hmax*hmax);
      #endif

      //perf_log.push("push back");
      stab_dofs_cols2.push_back(dof_indices_p[0]);
      stab_dofs_vals2.push_back(-factor);
      stab_dofs_cols2.push_back(neighbor_dof_indices_p[0]);
      stab_dofs_vals2.push_back(factor);
      //perf_log.pop("push back");
			}
		}

	//Lots of mallocs happen during the first assemble due to the unexpected sparsity pattern;
	//The implicit_neighbor_dof does not seem to work ?
	//perf_log.push("kstab");
  DenseMatrix<Number> Kstab2;
  Kstab2.resize(1, stab_dofs_vals2.size());
 	for (int i=0; i < stab_dofs_vals2.size(); i++) {
	 	Kstab2(0,i)=stab_dofs_vals2[i];
	}
 	std::vector<unsigned int> stab_dofs_rows2;
	stab_dofs_rows2.push_back(dof_indices_p[0]);
	//system.matrix->add_matrix(Kstab2,stab_dofs_rows2,stab_dofs_cols2);
  test(4);
 
	//perf_log.pop("kstab");
  #endif 
  //endif USE_STAB

#if USE_STAB 

  #if THREED
  #include "assemble_stokes_bcs_p1p1p0_anal_sine.cpp"
  #include "pin_pressure.cpp"
  #endif

  #if !THREED
  #include "assemble_stokes_bcs_p1p1p0_anal_sine.cpp"
  #include "square_pin_pressure.cpp"
  #endif


#endif

 // system.matrix->add_matrix (Ke, dof_indices);
  system.rhs->add_vector    (Fe, dof_indices);

} // end of element loop
  
    //system.matrix->close();
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


