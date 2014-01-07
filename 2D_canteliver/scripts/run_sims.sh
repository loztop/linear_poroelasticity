#!/bin/bash
# My first script

f_prefix="2dNEW_"


NE=(16 32)
NT=(8 16)

#clpc directory
exe_directory="/users/lorenzb/Dphil/libmesh_projetcs/linear_poro_2Dstabexp/"

#loztop directory
#exe_directory="/home/loztop/Dropbox/Dphil/libmesh_projetcs/linear_poro_2Dstabexp/"

matfiles_dir="data/matfiles/"
data_dir="data/"

res_directory_mat="$exe_directory$matfiles_dir"

res_directory_data="$exe_directory$data_dir"


exe_filename="ex11-opt"


str_nt='NT_'
str_ne='NE_'


for i in ${NE[@]}
do
for j in ${NT[@]}
do




output_file_name_mat="$res_directory_mat$f_prefix"_"$str_nt$j"_"$str_ne$i"_.mat" "

output_file_name_data="$res_directory_data$f_prefix"_"$str_nt$j"_"$str_ne$i"_" "

exe_str="$exe_directory$exe_filename $j $i $output_file_name_mat $output_file_name_data" 
   echo $exe_str

`$exe_str`

done
done



