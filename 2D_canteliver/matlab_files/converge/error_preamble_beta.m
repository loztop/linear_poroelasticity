%%Setup simulations
%f_prefix='2d_rescale_ft1_dtdivz_stab1dt'
clear all;

f_prefix='2dNEW_beta'

NE=[16,32];
NT=[8,16];
BETA=[-1 -0.5 0 0.5];

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

%%MPI run command, number of procs
mpirun=' mpirun -np 1 '

%First commpile all the sources
make_str_exe_error = strcat(['!make -j 1 all -C ' exe_directory]);

str_nt='NT_';
str_ne='NE_';
str_beta='BETA_';
