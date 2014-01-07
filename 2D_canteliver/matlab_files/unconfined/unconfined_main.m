clear all;
%close all;

tec='.tec';

%P2P2P1 solution
%base='/home/loztop/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cylinderP1k1_133sym_10Nt_';
%base='/home/loztop/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cylinderP1k1_133sym_20Nt_';


%base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_0p_133sym_40Nt_1T_';

%base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_0p_737sym_40Nt_20T_01k_';


base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_0p_133sym_20Nt_10T_01k_001stab_';

base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_np_133sym_20Nt_2T_005k_001stab_';

base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_np_133sym_20Nt_2T_005k_1stab_';

base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_np_133sym_20Nt_2T_005k_01stab_';

base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_np_133sym_20Nt_2T_005k_001stab_';

base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_np_133sym_10Nt_2T_005k_0001stab_';

base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_np_737sym_10Nt_2T_005k_0001stab_';

base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_np_5567sym_100Nt_10T_01k_0001stab_';

%base='/home/scratch/Dropbox/Dphil/libmesh_projetcs/nitsche_fix/data/cyl_classic_np_737sym_100Nt_10T_01k_0001stab_';



%Parameters for simulation (analytical solution)
params.T=5; %End Time
params.Nt=50;  %Number of time steps
params.dt=params.T/params.Nt;    %Size of time step

params.a=1;    %Radius of cylinder
params.v=0.15; %poisson ratio of elastic skeleton
params.E=1; %Youngs modulus of elastic skeleton

params.lambda=(params.E*params.v)/((1+params.v)*(1-2*params.v));   %elastic coefficent
params.mu=params.E/(2*(1+params.v));   %elastic coefficent

%params.Hk=1;   %aggrefate modulus of elastic skeleton (Hk=lambda + 2*mu !)
params.Hk=params.lambda+2*params.mu;   %aggrefate modulus of elastic skeleton (Hk=lambda + 2*mu !)

params.k=0.1;    %dynamic permeability
params.tg=1/(params.Hk*params.k/(params.a*params.a));  %characteristic time of diffusion
params.ez=0.05; %Amplitude of applied axial strain
params.Nodes=133; %Number of nodes

%create matrix to store results (NTx(number of idx x===1))
a=zeros(params.Nt,4);


for i=1:params.Nt
i
%Load time step data    
fname=strcat(base,num2str(i),tec);
importfile(fname);

header_info=textdata{3};
TEMP=data;
clear data;
data=TEMP(1:params.Nodes,:);

%This might change depending on what order vraibles are stored

%The reference positions
x=data(:,11);
y=data(:,12);
z=data(:,13);

%The displacements
s_u=data(:,14);
s_v=data(:,15);
s_w=data(:,16);

%Find index of nodes where x=1, y=1, z=1
idx_x_1=find(x==1);
idx_y_1=find(y==1);
idx_z_1=find(z==1);

%get the radial (x direction) displacemnet at points x=1
s_u_x_1=s_u(idx_x_1);

%Store these displacement for this time in a vector
%a(i,:)=s_u_x_1;

%Store the max and min values of radial displacement
%num_u(i)=max(s_u);
%num_u_min(i)=min(s_u);

num_u_min(i)=min(s_v);

num_u(i)=(max(s_v)-min(s_v))/2;


end
num_u

%Creat numeical x-axis and y solution vector
num_x=params.dt.*(1:1:params.Nt)./params.tg;
num_y=1.*(num_u-1)./(1*params.a*params.ez);

%Calculate analytical solution
[b_y,b_x]=bessel(params);

%plot numerical solution
figure;
%plot analytical solution
plot(b_x,b_y,'k','LineWidth',3);
hold all

plot(num_x,num_y,'rx','MarkerSize',14,'LineWidth',2);
hold all


save unconfined.mat
clear
load unconfined.mat

%axis([0 1.1 0 0.6])

% %Calculate Root-mean-square deviation of available co incideing time points
% rmsd=0;
% for i=1:params.Nt
% anal_y(i)=y(i*50);
% rmsd=rmsd+(anal_y(i)-num_y(i))^2;
% end
% rmsd=sqrt(rmsd/params.Nt)



