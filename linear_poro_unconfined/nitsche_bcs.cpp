for (unsigned int s=0; s<elem->n_sides(); s++)
{
   if (elem->neighbor(s) == NULL)
    {   
        AutoPtr<Elem> side (elem->build_side(s));






        //fluid face 
  			const std::vector<std::vector<Real> >&  phi_face_f =  fe_face_f->get_phi();
  			const std::vector<std::vector<RealGradient> >& dphi_face_f = fe_face_f->get_dphi();
  			const std::vector<Real>& JxW_face_f = fe_face_f->get_JxW();
  			const std::vector<Point>& qface_point_f = fe_face_f->get_xyz();
  			const std::vector<Point>& face_normals_f = fe_face_f->get_normals();
  			fe_face_f->reinit(elem,s);  

				 //pressure face 
				const std::vector<std::vector<Real> >&  phi_face_p =  fe_face_p->get_phi();
  			const std::vector<std::vector<RealGradient> >& dphi_face_p = fe_face_p->get_dphi();
  			const std::vector<Real>& JxW_face_p = fe_face_p->get_JxW();
  			const std::vector<Point>& qface_point_p = fe_face_p->get_xyz();
  			const std::vector<Point>& face_normals_p = fe_face_p->get_normals();
  			fe_face_p->reinit(elem,s); 

				

		// h elemet dimension to compute the interior penalty penalty parameter
        const unsigned int elem_b_order = static_cast<unsigned int> (fe_face_f->get_order());
        const double h_elem = elem->volume()/side->volume() * 1./pow(elem_b_order, 2.);

		Real hmax=(*elem).hmax();
		//std::cout<<" hmax "<<hmax<<std::endl;
		//std::cout<<" h_elem "<<h_elem<<std::endl;

		Real penalty =20;
		Real darcy_fac=1000*(1.0/(hmax*hmax));


    for (unsigned int qp=0; qp<qface_f->n_points(); qp++)
	{

		Real xf =	qface_point_f[qp](0);
        Real yf = qface_point_f[qp](1);
        #if THREED
        Real zf = qface_point_f[qp](2);
        #endif



	if( (qface_point_f[qp](0)<0.01) || (qface_point_f[qp](0)>0.99) ) {
		Real bc_value=4*(1-yf)*yf;
	//	std::cout<< "bc_value "<<bc_value<<std::endl;
     	bc_value=0;

			// (p,v.n)_{Gamma} term
          	for (unsigned int i=0; i<phi_face_p.size(); i++){
              	for (unsigned int j=0; j<phi_face_f.size(); j++){
								Kpx(i,j) += JxW_face_f[qp]*phi_face_p[i][qp]*phi_face_f[j][qp]*face_normals_f[qp](0);
								Kpy(i,j) += JxW_face_f[qp]*phi_face_p[i][qp]*phi_face_f[j][qp]*face_normals_f[qp](1);
								#if THREED
								Kpz(i,j) += JxW_face_f[qp]*phi_face_p[i][qp]*phi_face_f[j][qp]*face_normals_f[qp](2);
								#endif
		      	}
		    }


			// (penalty/h_elem*u.n,v.n)_{Gamma} term //stability term
			for (unsigned int i=0; i<phi_face_f.size(); i++)		{
              	for (unsigned int j=0; j<phi_face_f.size(); j++)			{
					Kxx(i,j) += penalty/h_elem*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](0)*phi_face_f[j][qp]*face_normals_f[qp](0);
					Kyy(i,j) += penalty/h_elem*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1)*phi_face_f[j][qp]*face_normals_f[qp](1);
					#if THREED
					Kzz(i,j) += penalty/h_elem*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](2)*phi_face_f[j][qp]*face_normals_f[qp](2);
					#endif
		    	}
		    }


            // RHS contributions
			
			//(u0.n,q)_D    //I think this is darcy ?!
			for (unsigned int i=0; i<phi_face_p.size(); i++){				
				Fp(i) +=  bc_value*JxW_face_f[qp]*face_normals_f[qp](0)*phi_face_p[i][qp];
				Fp(i) +=  0*JxW_face_f[qp]*face_normals_f[qp](1)*phi_face_p[i][qp];
				#if THREED
          		Fp(i) +=  0*JxW_face_f[qp]*face_normals_f[qp](2)*phi_face_p[i][qp];
				#endif
			}

			//(u0.n,penalty/h_elem*v.n)_D    
			for (unsigned int i=0; i<phi_face_p.size(); i++){				
				Fx(i) +=  bc_value*penalty/h_elem*JxW_face_f[qp]*face_normals_f[qp](0)*phi_face_f[i][qp]*face_normals_f[qp](0);
				Fy(i) +=  bc_value*penalty/h_elem*JxW_face_f[qp]*face_normals_f[qp](1)*phi_face_f[i][qp]*face_normals_f[qp](1);
				#if THREED
          		Fz(i) +=  bc_value*penalty/h_elem*JxW_face_f[qp]*face_normals_f[qp](2)*phi_face_f[i][qp]*face_normals_f[qp](2);
				#endif
			}

          
        }  //end if

	
	} //end qp
   } //if (elem->neighbor(s) == NULL)
}// end boundary condition section  


