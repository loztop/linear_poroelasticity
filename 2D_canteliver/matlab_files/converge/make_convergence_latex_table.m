clear all;
close all;

importfile('error_results.mat')

A=error_results;

%matrix = [1.5 1.764; 3.523 0.2];
%rowLabels = {'row 1', 'row 2'};
%columnLabels = {'col 1', 'col 2'};
%matrix2latex(matrix, 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
%matrix2latex(matrix, 'out.tex');
 


D1=zeros(8,8);
D2=zeros(8,8);

for j=1:length(A)    
    r=A(j,1)
    c=A(j,2)
    D1(r,c)= A(j,3); 
    D2(r,c)= A(j,4);   
end
    
D1(:,1)=[];
D1(3,:)=[];
D1([4:6],:)=[];

D2(:,1)=[];
D2(3,:)=[];
D2([4:6],:)=[];

%create row lables

NT=unique(A(:,1));
NE=unique(A(:,2));

for i=1:length(NT)
rowLabels(i) = {num2str(NT(i))};
end

for i=1:length(NE)
columnLabels(i) = {num2str(NE(i))};
end


matrix2latex(D1, 'conv_tab_L2.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'large');

matrix2latex(D2, 'conv_tab_H1.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'large');
