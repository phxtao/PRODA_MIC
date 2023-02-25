function fun_save_post_da_simu(save_pathway, save_var_namelist, save_var)
var_num = length(save_var_namelist);
for ivar = 1:var_num
    var_data_middle = save_var{ivar};
%     eval(['save_var_namelist{', num2str(ivar), '} = var_data_middle;']);
    
    save(save_pathway{ivar}, 'var_data_middle');
end

end
