function all_practCE = get_practCE(micro_variables, macro_variables, tau, method);
	
	% gets as inputs two structures - one with micro, and one with macro variables - 
	% where each cell is a micro or macro variable;
	% returns a structure where cells denote different combinations of micro and macro variables; 
	% each cell, in turn, is a structure with practical CE, DC, & CD
	
	% get number of micro and macro variables, respectively
	n_micro_variables = length(fieldnames(micro_variables)); 
	n_macro_variables = length(fieldnames(macro_variables)); 
	
	% get names of micro and macro variables, respectively
	fieldnames_macro = fieldnames(macro_variables);
	fieldnames_micro = fieldnames(micro_variables);
	
	% loop over all micro and macro variables and calculate practical CE
	for i = 1:n_micro_variables;
		micro = micro_variables.(fieldnames_micro{i});
		
		for j = 1:n_macro_variables;
			macro = macro_variables.(fieldnames_macro{j});

			% calculate practical CE, DC, and CE for a given micro-macro combination, and store 
			% results in struct
			practCE.pract_ce = EmergencePsi(micro', macro', tau, method);
			practCE.pract_dc = EmergenceDelta(micro', macro', tau, method);
			practCE.pract_cd = EmergenceGamma(micro', macro', tau, method);
			
			% store struct in another struct with micro-macro combination as fieldname
			all_practCE.([fieldnames_micro{i} '_' fieldnames_macro{j}]) = practCE;
		end 
	end 
	
end 