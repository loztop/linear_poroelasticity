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
            
//block flow bottom in z direction 
#if THREED
          if ((elem->node(n) == side->node(ns)) && ( (zf<0.001) ))
#endif           
          {
            int source_dof = node->dof_number(system.number(), z_var, 0);
            

            if((source_dof<12345678) && (source_dof>-1)){
             rows_values.push_back(0);
             rows.push_back(source_dof);
            }


            
          }  //end if



          //block the top
          #if THREED
          if ((elem->node(n) == side->node(ns)) && (  (zf>0.999) ))
          #endif
          {
          
          
            int source_dof = node->dof_number(system.number(), z_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(0);
              rows.push_back(source_dof);
            }
            }  //end if


           //Constrain pressure to zero at boundary
        
          if ((elem->node(n) == side->node(ns)) && ( sqrt(xf*xf + yf*yf)  >0.99)  )
          {

            int source_dof = elem->dof_number(system.number(), p_var, 0);

           // std::cout<< (*node) <<std::endl;
           // std::cout<< source_dof <<std::endl;

            if((source_dof<12345678) && (source_dof>-1)){
           //   pressure_rows_values.push_back(0);
           //   pressure_rows.push_back(source_dof);
            }
          }  //end if



          }
        }
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


