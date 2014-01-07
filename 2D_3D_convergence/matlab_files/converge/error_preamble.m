%%Setup simulations
f_prefix='2D_stab_0.95'

NE=[8, 16, 32];
NT=[8,16, 32];

%%Set up the correct directories

%loztop
exe_directory = '/home/loztop/Dropbox/Dphil/libmesh_projetcs/poro_paper_sims/2D_convergence_delta2/'
%clpc59
%exe_directory = '/users/lorenzb/Dphil/libmesh_projetcs/poro_paper_sims/2D_convergence_delta2/'

exe_filename = 'ex11-opt'

res_directory =strcat([exe_directory 'data/matfiles/'])
res_directory_plot =strcat([exe_directory 'data/plots/'])

%at home
%res_directory = '/home/loztop/Dropbox/Dphil/libmesh_projetcs/linear_poro_convergence_2Dstabexp/data/matfiles/'
%res_directory_plot = '/home/loztop/Dropbox/Dphil/libmesh_projetcs/linear_poro_convergence_2Dstabexp/data/plots/'

%%Run simulations
%run from commandline
% matlab -nojvm -nodisplay -nosplash -r "run_sims"
%export LD_PRELOAD=/usr/lib64/libstdc++.so.6 matlab -nojvm -nodisplay -nosplash -r "run_sims"

%%MPI run command, number of procs
mpirun=' mpirun -np 1 '

%First commpile all the sources
make_str_exe_error = strcat(['!make -j 1 all -C ' exe_directory]);


str_nt='NT_';
str_ne='NE_';

fine_NT=NT(end);
fine_NE=NE(end);

fine_file_name = strcat([res_directory f_prefix,...
     '_' str_nt int2str(fine_NT) '_' str_ne int2str(fine_NE) '_.xdr']);
 