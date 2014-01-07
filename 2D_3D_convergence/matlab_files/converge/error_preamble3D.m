%%Setup simulations
%f_prefix='stokes_2D'
f_prefix='3Dasymtimeloz0_001stab'

NE=[4,8,12];
NT=[4,8,12];

%%Set up the correct directories


%loztop
exe_directory = '/home/loztop/Dropbox/Dphil/libmesh_projetcs/linear_poro_convergence_june25/'

exe_filename = 'ex11-opt'

%clpc59
res_directory =strcat([exe_directory 'data/matfiles/'])
res_directory_plot =strcat([exe_directory 'data/plots/'])
%%Run simulations
%run from commandline
% matlab -nojvm -nodisplay -nosplash -r "run_sims"
%export LD_PRELOAD=/usr/lib64/libstdc++.so.6 matlab -nojvm -nodisplay -nosplash -r "run_sims"

%%MPI run command, number of procs
mpirun=' mpirun -np 1 '


str_nt='NT_';
str_ne='NE_';

fine_NT=NT(end);
fine_NE=NE(end);

fine_file_name = strcat([res_directory f_prefix,...
     '_' str_nt int2str(fine_NT) '_' str_ne int2str(fine_NE) '_.xdr']);
 