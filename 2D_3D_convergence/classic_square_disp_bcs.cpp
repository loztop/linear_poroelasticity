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
 
/*           
//Fix displacement at x==0
          if ((elem->node(n) == side->node(ns)) && ( (xf<0.001) ))
          {
            int source_dof = node->dof_number(system.number(), u_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(0);
              rows.push_back(source_dof);
            }

            source_dof = node->dof_number(system.number(), v_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
             rows_values.push_back(0);
             rows.push_back(source_dof);
            }

            #if THREED
            source_dof = node->dof_number(system.number(), w_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
             rows_values.push_back(0);
             rows.push_back(source_dof);
            }
            #endif
            
          }  //end if
*/


//Prescribe displacement at x==0
          if ((elem->node(n) == side->node(ns)) &&  (xf<0.01))
          {
            int source_dof = node->dof_number(system.number(), u_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(0);
              rows.push_back(source_dof);
            }

          }  //end if


				//Prescribe displacement at x==1
          if ((elem->node(n) == side->node(ns)) &&  (xf>0.999))
          {
            int source_dof = node->dof_number(system.number(), u_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(-0.4*progress);
              rows.push_back(source_dof);
            }

          }  //end if



          }
        }
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


