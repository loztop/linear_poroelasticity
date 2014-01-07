  const Real dt    = es.parameters.get<Real>("dt");
  const Real DELTA    = es.parameters.get<Real>("DELTA");

 // const Real time    = es.parameters.get<Real>("time");
 // const Real progress    = es.parameters.get<Real>("progress");
 // const Real delta_power    = es.parameters.get<Real>("beta");


//Real mu = E/(2*(1+NU));
//Real lambda = (E*NU)/((1+NU)*(1-2*NU));
//Real elastic_fac=pow(dt,0);
//Real fluid_gradp_fac=pow(dt,0);


  PerfLog perf_log("Assemble");
  libmesh_assert (system_name == "Last_non_linear_soln");
 
  const MeshBase& mesh = es.get_mesh();
  
  const unsigned int dim = mesh.mesh_dimension();
  
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem> ("Last_non_linear_soln");

  const unsigned int u_var = system.variable_number ("s_u");
  const unsigned int v_var = system.variable_number ("s_v");
  #if THREED
  const unsigned int w_var = system.variable_number ("s_w");
  #endif
  const unsigned int p_var = system.variable_number ("s_p");
  const unsigned int x_var = system.variable_number ("x");
  const unsigned int y_var = system.variable_number ("y");
   #if THREED
  const unsigned int z_var = system.variable_number ("z");
  #endif
  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_disp_type = system.variable_type(u_var);
  FEType fe_vel_type = system.variable_type(x_var);

  // Get the Finite Element type for "p".
  FEType fe_pres_type = system.variable_type(p_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  AutoPtr<FEBase> fe_disp  (FEBase::build(dim, fe_disp_type));
  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
    
  // Build a Finite Element object of the specified type for
  // the pressure variables.
  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  
  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_disp->attach_quadrature_rule (&qrule);
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  
  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe_vel->get_JxW();
  
  const std::vector<Point>& q_point = fe_vel->get_xyz();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi = fe_disp->get_dphi();
  const std::vector<std::vector<Real> >& phi = fe_disp->get_phi();
  const std::vector<std::vector<RealGradient> >& f_dphi = fe_vel->get_dphi();
  const std::vector<std::vector<Real> >& f_phi = fe_vel->get_phi();

  // The element shape functions for the pressure variable
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();


#if NITSCHE
	AutoPtr<FEBase> fe_face_d (FEBase::build(dim, fe_disp_type));          
	AutoPtr<QBase> qface_d(fe_disp_type.default_quadrature_rule(dim-1));
	fe_face_d->attach_quadrature_rule (qface_d.get());

	AutoPtr<FEBase> fe_face_f (FEBase::build(dim, fe_vel_type));          
	AutoPtr<QBase> qface_f(fe_vel_type.default_quadrature_rule(dim-1));
	fe_face_f->attach_quadrature_rule (qface_f.get());

	AutoPtr<FEBase> fe_face_p (FEBase::build(dim, fe_pres_type));          
	AutoPtr<QBase> qface_p(fe_pres_type.default_quadrature_rule(dim-1));
	fe_face_p->attach_quadrature_rule (qface_p.get());

	//AutoPtr<FEBase> fe_face_ref (FEBase::build(dim, fe_vel_type_ref));          
	//AutoPtr<QBase> qface_ref(fe_vel_type_ref.default_quadrature_rule(dim-1));
	//fe_face_ref->attach_quadrature_rule (qface_ref.get());
#endif

  const DofMap & dof_map = system.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseMatrix<Number> Kstab;

#if !THREED
  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke), Kux(Ke), Kuy(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke), Kvx(Ke), Kvy(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke), Kpx(Ke), Kpy(Ke),
    Kxu(Ke), Kxv(Ke), Kxp(Ke), Kxx(Ke), Kxy(Ke),
    Kyu(Ke), Kyv(Ke), Kyp(Ke), Kyx(Ke), Kyy(Ke);
#endif

#if THREED
  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke), Kux(Ke), Kuy(Ke), Kuz(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke), Kvx(Ke), Kvy(Ke), Kvz(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke), Kwx(Ke), Kwy(Ke), Kwz(Ke),
    Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke), Kpx(Ke), Kpy(Ke), Kpz(Ke),
    Kxu(Ke), Kxv(Ke), Kxw(Ke), Kxp(Ke), Kxx(Ke), Kxy(Ke), Kxz(Ke),
    Kyu(Ke), Kyv(Ke), Kyw(Ke), Kyp(Ke), Kyx(Ke), Kyy(Ke), Kyz(Ke),
    Kzu(Ke), Kzv(Ke), Kzw(Ke), Kzp(Ke), Kzx(Ke), Kzy(Ke), Kzz(Ke);
#endif

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe),
    Fx(Fe),
    Fy(Fe);

  #if THREED
  DenseSubVector<Number>
    Fw(Fe),
    Fz(Fe);
  #endif

  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_p;
  std::vector<unsigned int> dof_indices_x;
  std::vector<unsigned int> dof_indices_y;
  #if THREED
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_z;
  #endif
  
 MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  
std::vector<  int > rows;
std::vector<  int > pressure_rows;
std::vector< Real > rows_values;
std::vector< Real > pressure_rows_values;

std::vector<unsigned int> stab_dofs_rows;
std::vector<unsigned int> stab_dofs_cols;
std::vector<Real> stab_dofs_vals;
