clear all;

load unconfined.mat

hFig=figure;
set(hFig, 'Position', [100 100 900 500])
set([gca]             , ...
    'FontSize'   , 12           );



%plot analytical solution
%plot(b_x,b_y,'k','LineWidth',3);
%hold all

%plot(num_x,num_y,'rx','MarkerSize',14,'LineWidth',2);
%hold all



num_hand=plot(num_x,num_y,'-ko','LineWidth',2,'MarkerSize',7)
hold on;
an_hand=plot(b_x,b_y,'-k','LineWidth',2,'MarkerSize',8)
hold on;


%Add the legend and labels
title('Unconfined compression relaxation test','interpreter','latex','FontSize',15);
xlabel('Nondimensional time $(t/t_{g})$ ','interpreter','latex','FontSize',19)
ylabel('Radial displacement $(u/a\epsilon_{0})$ ','interpreter','latex','FontSize',19)

hLegend = legend( ...
   [num_hand,an_hand], ...
  'Numerical solution' ,'Analytical solution',...
  'location', 'NorthEast' );

h=hFig;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

output_plot_filename='~/Dropbox/Dphil/linear_poro_paper/diagrams/unconfined_results'
print(h,output_plot_filename,'-dpdf','-r0')

