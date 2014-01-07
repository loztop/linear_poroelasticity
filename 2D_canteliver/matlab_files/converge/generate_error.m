%%include the preamble file setup stuff
error_preamble;

%eval(make_str_exe_error) - have to make on command line.
for(j=1:length(NE))

for(i=1:length(NT))

 output_file_name = strcat([res_directory f_prefix,...
     '_' str_nt int2str(NT(i)) '_' str_ne int2str(NE(j)) '_.mat'] );
  
exec_str = strcat(['!' exe_directory exe_filename,...
     ' ' int2str(NT(i)) ' ' int2str(NE(j)),...
    ' ' output_file_name])

eval(exec_str) 

end
end




