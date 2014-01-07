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

    for (unsigned int qp=0; qp<qface_f->n_points(); qp++)
	{

		Real xf =	qface_point_f[qp](0);
        Real yf = qface_point_f[qp](1);
        #if THREED
        Real zf = qface_point_f[qp](2);
        #endif

            // RHS contributions
			if( (yf>0.99) ){

			Point disp_value;
			disp_value(0)=0;
			disp_value(1)=0;
			disp_value(2)=0;

				//-(u0.n,dnv)_D    //consistency
				for (unsigned int i=0; i<phi_face_f.size(); i++){				
					
					Fv(i) -= 1*JxW_face_f[qp] *phi_face_f[i][qp]; 
			
				}			
		}
         	
	} //end qp
   } //if (elem->neighbor(s) == NULL)
}// end boundary condition section  


