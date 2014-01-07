clear all;

tec='.tec';
base='/auto/users/lorenzb/Dphil/libmesh_projetcs/linear_poro_canteliver/data/cant_nopin_disc_96_';
%base='/auto/users/lorenzb/Dphil/libmesh_projetcs/linear_poro_canteliver/data/cant_nopin_disc_unstable_';

fname=strcat(base,num2str(5),tec);

importfile(fname);
textdata{3}
x=data(:,9);
y=data(:,10);

s_p=data(:,6);

idx_x_1=find(x==0.25);
idx_x_2=find(x==0.5);
idx_x_3=find(x==0.75);

s_p_x_1=s_p(idx_x_1);
s_p_x_2=s_p(idx_x_2);
s_p_x_3=s_p(idx_x_3);

x_axis=[1:size(s_p_x_1,1)-1]./(size(s_p_x_1,1)-1)
hFig=figure;
h_hand_1=plot(x_axis,s_p_x_1(1:end-1),'-k','LineWidth',2,'MarkerSize',8)
hold on;
h_hand_2=plot(x_axis,s_p_x_2(1:end-1),'--k','LineWidth',2,'MarkerSize',8)
hold on;
h_hand_3=plot(x_axis,s_p_x_3(1:end-1),'-xk','LineWidth',1,'MarkerSize',8)

%Add the legend and labels
% title('Cantilever bracket problem','interpreter','latex','FontSize',17);
% 
 xlabel('$y$','interpreter','latex','FontSize',17)
 ylabel('Pressure','interpreter','latex','FontSize',17)

 hLegend = legend( ...
  [h_hand_1,h_hand_2,h_hand_3], ...
  strcat(['x=0.25 ' ]) , ...
  strcat(['x=0.5 ' ]), ...
  strcat(['x=0.75 ']) ,...
  'location', 'SouthWest' );

h=hFig;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%output_plot_filename='~/Dropbox/Dphil/linear_poro_paper/diagrams/convergence_rescale_2d_asym'
output_plot_filename='/auto/users/lorenzb/Dphil/linear_poro_paper/diagrams/cantilever_plot';

print(h,output_plot_filename,'-dpdf','-r0')
