%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% change_fieldnames_in_table.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a quick-and-dirty script that changes the fieldname 'coupling_mag' to 'global_couplingmag'
% in each struct in the table. This is only possible by building a new
% table and reassigning it to the original table name.
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_cross_coup_mag = 2.^linspace(cross_coup_range(1),cross_coup_range(2), n_cross_coup); 
all_rmi   = 2.^linspace(rmi_range(1), rmi_range(2),  n_rmi  ); 

% row & column values for table
all_cross_coup_mag_str = {};
for t = 1:length(all_cross_coup_mag)
        all_cross_coup_mag_str{t} = num2str(all_cross_coup_mag(t));
end

all_rmi_str = {};
for e = 1:length(all_rmi)
        all_rmi_str{e} = num2str(all_rmi(e));
end

for i=1:size(results_table,1)
    for j=1:size(results_table,2)
        current_struct = table2array(results_table(i,j));
        fnames = fieldnames(current_struct);
        old_idx = strcmp(fnames, 'coupling_mag');
        fnames{old_idx} = 'global_coupling_mag';
        current_struct = cell2struct(struct2cell(current_struct), fnames);
        new_struct(i,j) = current_struct;
    end
end

if strcmp(measure_in_filename, 'phiid_measures_mmi') || strcmp(measure_in_filename, 'phiid_measures_ccs')
	part_of_filename = measures_in_table;
else 
	part_of_filename = measure_in_filename;
end

% as table
results_table = array2table(new_struct, 'RowNames', all_cross_coup_mag_str, 'VariableNames', all_rmi_str);
results_table_filepath = fullfile(pathout_data_measures, [model_specific_filename '_' part_of_filename '.mat']);
fprintf('\nSaving results table file ''%s'' ...', results_table_filepath);
save(results_table_filepath, 'results_table', '-v7.3');
fprintf(' d1.0e+07 one\n');

clear new_struct;
clear results_table;

%save('/media/nadinespy/NewVolume1/work/phd/projects/mec_experiments/mec_simulations/results/analyses/multivar_autoreg/rand_cross_coup_rmi/8mvar_lag1/8mvar_lag1_rand_cross_coup_rmi_phiid_mmi_uc_table.mat', 'results_table', '-v7.3');
