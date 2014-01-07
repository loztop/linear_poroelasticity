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

void assemble_stiffness (EquationSystems& es,
                      const std::string& system_name)
{
#include "assemble_preamble.cpp"

std::cout<< "Assemble Stiffness "<<std::endl;

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
       
        for (unsigned int i=0; i<n_u_dofs; i++){
          for (unsigned int j=0; j<n_u_dofs; j++){
            Kuu(i,j) += JxW[qp]*dphi[i][qp]*dphi[j][qp];
            Kvv(i,j) += JxW[qp]*dphi[i][qp]*dphi[j][qp];
          }
        }


          // up coupling
          for (unsigned int i=0; i<n_u_dofs; i++){
            for (unsigned int j=0; j<n_p_dofs; j++){
              Kup(i,j) += -JxW[qp]*dphi[i][qp](0)*psi[j][qp];
              Kvp(i,j) += -JxW[qp]*dphi[i][qp](1)*psi[j][qp];
            }
          }


          Real fluid_fac=(1.0/KPERM);

          //Mass Matrix
          for (unsigned int i=0; i<n_x_dofs; i++)
            for (unsigned int j=0; j<n_x_dofs; j++)
              Kxx(i,j) += fluid_fac*JxW[qp]*(f_phi[i][qp]*f_phi[j][qp]);

          for (unsigned int i=0; i<n_y_dofs; i++)
            for (unsigned int j=0; j<n_y_dofs; j++)
              Kyy(i,j) += fluid_fac*JxW[qp]*(f_phi[i][qp]*f_phi[j][qp]);


          //Grad P
          for (unsigned int i=0; i<n_x_dofs; i++){
            for (unsigned int j=0; j<n_p_dofs; j++){
              Kxp(i,j) += -JxW[qp]*f_dphi[i][qp](0)*psi[j][qp];
              Kyp(i,j) += -JxW[qp]*f_dphi[i][qp](1)*psi[j][qp];
            }
          }

            //Elastcity (Mixture) momentum equation

#if FULL_ELASTIC 
        Real fac=2*elastic_fac;
          for (unsigned int i=0; i<n_u_dofs; i++){
            for (unsigned int j=0; j<n_u_dofs; j++){
              Kuu(i,j) += fac*mu*JxW[qp]*(dphi[i][qp](0)*dphi[j][qp](0));
              Kuu(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](1)*dphi[j][qp](1));
              Kuu(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](2)*dphi[j][qp](2));
              Kuv(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](1)*dphi[j][qp](0));
              #if THREED
              Kuw(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](2)*dphi[j][qp](0));
              #endif

              Kvv(i,j) += fac*mu*JxW[qp]*(dphi[i][qp](1)*dphi[j][qp](1));
              Kvv(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](0)*dphi[j][qp](0));
              Kvv(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](2)*dphi[j][qp](2));
              Kvu(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](0)*dphi[j][qp](1));
              #if THREED
              Kvw(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](2)*dphi[j][qp](1));
              #endif

              #if THREED
              Kww(i,j) += fac*mu*JxW[qp]*(dphi[i][qp](2)*dphi[j][qp](2));
              Kww(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](0)*dphi[j][qp](0));
              Kww(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](1)*dphi[j][qp](1));
              Kwu(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](0)*dphi[j][qp](2));
              Kwv(i,j) += fac*mu*JxW[qp]*(1.0/2.0)*(dphi[i][qp](1)*dphi[j][qp](2));
              #endif
        }
    }
            //Divergence term
          Real lam_fac=elastic_fac;
          elastic_fac=1;
          for (unsigned int i=0; i<n_u_dofs; i++){
            for (unsigned int j=0; j<n_u_dofs; j++){
              Kuu(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](0)*dphi[j][qp](0);
              Kuv(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](0)*dphi[j][qp](1);
              #if THREED
              Kuw(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](0)*dphi[j][qp](2);
              #endif

              Kvu(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](1)*dphi[j][qp](0);
              Kvv(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](1)*dphi[j][qp](1);
              #if THREED
              Kvw(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](1)*dphi[j][qp](2);
              #endif

              #if THREED
              Kwu(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](2)*dphi[j][qp](0);
              Kwv(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](2)*dphi[j][qp](1);
              Kww(i,j) += lam_fac*lambda*JxW[qp]*dphi[i][qp](2)*dphi[j][qp](2);
              #endif
            }
          }
#endif


          //Mass conservation of mixture
          for (unsigned int i=0; i<n_p_dofs; i++){
            for (unsigned int j=0; j<n_u_dofs; j++){
              Kpu(i,j) +=  JxW[qp]*psi[i][qp]*dphi[j][qp](0);
              Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
             #if THREED
              Kpw(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](2);
              #endif
            }
          }

          for (unsigned int i=0; i<n_p_dofs; i++){
            for (unsigned int j=0; j<n_x_dofs; j++){
              Kpx(i,j) += dt*JxW[qp]*psi[i][qp]*f_dphi[j][qp](0);
              Kpy(i,j) += dt*JxW[qp]*psi[i][qp]*f_dphi[j][qp](1);
          #if THREED
              Kpz(i,j) += dt*JxW[qp]*psi[i][qp]*f_dphi[j][qp](2);
          #endif
            }
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

      //perf_log.push("push back");
      stab_dofs_cols2.push_back(dof_indices_p[0]);
      stab_dofs_vals2.push_back(-factor);
      stab_dofs_cols2.push_back(neighbor_dof_indices_p[0]);
      stab_dofs_vals2.push_back(factor);
			}
		}

	//Add stabilisation contr
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
	system.matrix->add_matrix(Kstab2,stab_dofs_rows2,stab_dofs_cols2);
	//perf_log.pop("kstab");

  #endif 
  //endif USE_STAB


  system.matrix->add_matrix (Ke, dof_indices);

} // end of element loop
  
    system.matrix->close();

  return;
}


