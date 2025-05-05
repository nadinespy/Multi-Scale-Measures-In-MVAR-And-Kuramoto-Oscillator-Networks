%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% change_fieldnames_in_table.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a quick-and-dirty script that changes the fieldname 'coupling_mag' to 
% 'global_couplingmag' in each struct in the table. This is only possible by building 
% a new table and reassigning it to the original table name. The new structure is 
% saved as a structure (and not as a table).
%
% ADD MORE DOCUMENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

old_fieldname = 'global_coupling_mag';
new_fieldname = 'global_coup';

for i=1:size(results_table,1)
    for j=1:size(results_table,2)
        current_struct = table2array(results_table(i,j));
        fnames = fieldnames(current_struct);
        old_idx = strcmp(fnames, old_fieldname);
        fnames{old_idx} = new_fieldname;
        current_struct = cell2struct(struct2cell(current_struct), fnames);
        new_struct(i,j) = current_struct;
    end
end

if strcmp(measure_in_filename, 'phiid_measures_mmi') || strcmp(measure_in_filename, ...
	'phiid_measures_ccs')
	part_of_filename = measures_in_table;
else 
	part_of_filename = measure_in_filename;
end

% save
results = new_struct;
clear new_struct;

results_filepath = fullfile(pathout_data_measures, [model_specific_filename ...
	'_' part_of_filename{1} '.mat']);
fprintf('\nSaving results file ''%s'' ...', results_filepath);
save(results_filepath, 'results', '-v7.3');

clear new_struct;
%clear results;
