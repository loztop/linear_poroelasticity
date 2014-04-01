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

          if ((elem->node(n) == side->node(ns)) && (  (xf>0.999) ))
          {
            int source_dof = node->dof_number(system.number(), u_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(-0.4*progress);
              rows.push_back(source_dof);
            }

            }  //end if


        //No outflow on constrained sides (x==0, x==1)
          if ((elem->node(n) == side->node(ns)) && ( (xf>0.999) || (xf<0.001)))
          {
            int source_dof = node->dof_number(system.number(), x_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(0);
              rows.push_back(source_dof);
            }
          }  //end if
          


      //p=0 dirichlet bc on outflow boundaries
          #if THREED
          if ((elem->node(n) == side->node(ns)) && (  (zf>0.99) ||zf<(0.01) || yf>(0.99) ||yf<(0.01)   ))
          #endif
          #if !THREED
          if ((elem->node(n) == side->node(ns)) && ( yf>(0.99) ||yf<(0.01)   ))
          #endif
          {

          #if !PRES_STAB
            int source_dof = node->dof_number(system.number(), p_var, 0);
          #endif

          #if PRES_STAB
            int source_dof = elem->dof_number(system.number(), p_var, 0);
          #endif
            //std::cout<< (*node) <<std::endl;
           // std::cout<< source_dof <<std::endl;
            if((source_dof<12345678) && (source_dof>-1)){
              pressure_rows_values.push_back(0);
              pressure_rows.push_back(source_dof);
            }
          }  //end if


          }
        }
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


