clear all;
close all;

error_preamble;
%NE=[NE,NE(end)];

 for(j=1:length(NE))

error_file_name = strcat([res_directory_plot f_prefix,...
     '_' str_nt int2str(NT(end)) '_' str_ne int2str(NE(j)) '_.mat']);

A = load('-ascii', error_file_name)


L2_u_h(j,1:2)=A(1:2);  
L2_u_h(j,3)=A(3);
H1_u_h(j,1:2)=A(1:2); 
H1_u_h(j,3)=A(4);

L2_w_h(j,1:2)=A(1:2);  
L2_w_h(j,3)=A(5);
Hd_w_h(j,1:2)=A(1:2);  
Hd_w_h(j,3)=A(6);

L2_p_h(j,1:2)=A(1:2);  
L2_p_h(j,3)=A(7);

 end
 
 %%
%%Do actual plotting
%close all;
%setup figure

hFig=figure;
set(hFig, 'Position', [100 100 600 600])
set([gca]             , ...
    'FontSize'   , 12           );

x_axis_NE=log(NE(1:end));

L2_u_h_hand=plot(x_axis_NE,log(L2_u_h(1:end,3)),'-kx','LineWidth',2,'MarkerSize',8)
hold on;
H1_u_h_hand=plot(x_axis_NE,log(H1_u_h(1:end,3)),'--kx','LineWidth',2,'MarkerSize',8)
hold on;

L2_w_h_hand=plot(x_axis_NE,log(L2_w_h(1:end,3)),'-b^','LineWidth',2,'MarkerSize',8);
hold on;
Hd_w_h_hand=plot(x_axis_NE,log(Hd_w_h(1:end,3)),'--b^','LineWidth',2,'MarkerSize',8)
hold on;


L2_p_h_hand=plot(x_axis_NE,log(L2_p_h(1:end,3)),'-ro','LineWidth',2,'MarkerSize',8);
hold on;

% %Add triangles
% polyfit(x_axis_NE,log(H1_h(1:end-1,3))',1)
% %Plot convergenc rate triangle
% %L2
% plot([1, 1.5, 1.5, 1],[-3,-3,-4.45,-3],'k-');
% hold on;
% text(1.55,-3.5,'2.9');
% hold on;
% 
 %H1
% plot([1, 1.5, 1.5, 1],[0,0,-0.8,0],'k-');
% hold on;
 %text(1.55,-0.5,'1.6');

%Add the legend and labels
title('2D h convergence','interpreter','latex','FontSize',19);
xlabel('$\log(1/h)$ ','interpreter','latex','FontSize',19)
ylabel('$\log(\mbox{error})$ ','interpreter','latex','FontSize',19)
grid on;
 set(gca,'YTick',[-8:1:4])
 set(gca,'XTick',[1:1:5])

hLegend = legend( ...
   [L2_u_h_hand,H1_u_h_hand,L2_w_h_hand,Hd_w_h_hand,L2_p_h_hand], ...
  'L2 displacement error' ,'H1 displacement error','L2 fluid velocity error' ,'Hdiv fluid velocity error','L2 pressure error',...
  'location', 'SouthWest' );

h=hFig;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

output_plot_filename='~/Dropbox/Dphil/linear_poro_paper/diagrams/convergence_2Dh_2'
print(h,output_plot_filename,'-dpdf','-r0')


