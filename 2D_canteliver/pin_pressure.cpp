    for (unsigned int s=0; s<elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {   
        AutoPtr<Elem> side (elem->build_side(s));

        for (unsigned int ns=0; ns<side->n_nodes(); ns++)
        {

          for (unsigned int n=0; n<elem->n_nodes(); n++)
          {
            Node *node = elem->get_node(n);
            const Real xf = (*node)(0);
            Real yf = (*node)(1);
            #if THREED
            Real zf = (*node)(2);
            #endif
            
            Number value_u = exact_2D_solution_u((*node), es.parameters,"null","void");
            Number value_v = exact_2D_solution_v((*node), es.parameters,"null","void");
						Number value_x = exact_2D_solution_x((*node), es.parameters,"null","void");
            Number value_y = exact_2D_solution_y((*node), es.parameters,"null","void");
            Number value_p = exact_2D_solution_p((*node), es.parameters,"null","void");
            #if THREED
            Number value_w = exact_2D_solution_w((*node), es.parameters,"null","void");
            Number value_z = exact_2D_solution_z((*node), es.parameters,"null","void");
            #endif

        
          #if !THREED
          if ((elem->node(n) == side->node(ns)) && ( (xf<0.001) && (yf<0.001) ))
          #endif 
          #if THREED
          if ((elem->node(n) == side->node(ns)) && ( (xf<0.001) && (yf<0.001) && (zf<0.001)  ))
          #endif 
          {
            int source_dof = elem->dof_number(system.number(), p_var, 0);

					//std::cout<<"pressure dof "<< source_dof <<" xf "<<xf<<" yf "<<yf <<std::endl;

            if((source_dof<12345678) && (source_dof>-1) ){
                pressure_rows_values.push_back(value_p);
                pressure_rows.push_back(source_dof);
            }
            
          }  //end if

          }
        }
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


