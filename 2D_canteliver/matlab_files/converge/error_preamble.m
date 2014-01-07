%%Setup simulations
%f_prefix='2d_rescale_ft1_dtdivz_stab1dt'
f_prefix='2dNEW'

%NT=[8,16,32,64,128];
%NE=[8,16,32,64,128];

%NT=[8,12,16,20,32,40];
%NE=[8,12,16,20,32,40];
%NT=ones(1,size(NT,2))*128;

NE=[16,32];
NT=[8,16];

S=1;


%%Set up the correct directories

%clpc
exe_directory = '/users/lorenzb/Dphil/libmesh_projetcs/linear_poro_2Dstabexp/'

%loztop
%exe_directory = '/home/loztop/Dropbox/Dphil/libmesh_projetcs/linear_poro_2Dstabexp/'

exe_filename = 'ex11-opt'

%clpc59
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
 