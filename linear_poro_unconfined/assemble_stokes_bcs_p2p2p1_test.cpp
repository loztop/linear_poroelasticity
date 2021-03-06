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

          if ((elem->node(n) == side->node(ns)) && ((yf>0.999)) ) {
            int source_dof = node->dof_number(system.number(), 0, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(0);
              rows.push_back(source_dof);
            }
            source_dof = node->dof_number(system.number(), 1, 0);
            if((source_dof<12345678) && (source_dof>-1)){
             rows_values.push_back(0);
             rows.push_back(source_dof);
            }
          }  //end if

          if ((elem->node(n) == side->node(ns)) && ((yf<-0.999)) ) {
            int source_dof = node->dof_number(system.number(), 1, 0);
            if((source_dof<12345678) && (source_dof>-1)){
             rows_values.push_back(0.3);
             rows.push_back(source_dof);
            }
          }  //end if

/*
      		if ((elem->node(n) == side->node(ns)) && ( (xf<-0.999) || (xf>0.999) || (yf>0.999) || (yf<-0.999) ))
          {
            
            int source_dof = node->dof_number(system.number(), 3, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(value_x);
              rows.push_back(source_dof);
            }

             source_dof = node->dof_number(system.number(), 4, 0);
            if((source_dof<12345678) && (source_dof>-1)){
              rows_values.push_back(value_y);
              rows.push_back(source_dof);
            }
          }  //end if
          */


/*
          if ((elem->node(n) == side->node(ns)) && ( (xf<-0.999) || (xf>0.999) || (yf>0.999) || (yf<-0.999) ))
          {
            int source_dof = elem->dof_number(system.number(), 2, 0);
            if((source_dof<12345678) && (source_dof>-1) ){
                pressure_rows_values.push_back(value_p);
                pressure_rows.push_back(source_dof);
            }
            
          }  //end if
          */

          }
        }
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


