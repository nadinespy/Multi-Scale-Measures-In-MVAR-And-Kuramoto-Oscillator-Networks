function all_DD = get_DD(micro_variables, macro_variables, method, tau_history, tau_steps, kraskov_param)

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
	
	if method == 'Kraskov'
		teCalc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorMultiVariateKraskov');
		teCalc.setProperty('k', num2str(kraskov_param));
	elseif method == 'Gaussian';
		teCalc = javaObject('infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorMultiVariateGaussian');
	elseif method == 'Discrete';
		teCalc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', 4, tau_history);
		% 4: base - number of symbols for each variable. E.g. binary variables are in base-2.
	end 

	% loop over all micro and macro variables and calculate practical CE
	for i = 1:n_micro_variables;
		micro = micro_variables.(fieldnames_micro{i});
		
		for j = 1:n_macro_variables;
			macro = macro_variables.(fieldnames_macro{j});
	
			if strcmp(method, 'Kraskov') | strcmp(method, 'Gaussian');
				micro_dim = size(micro, 2);
				macro_dim = size(macro, 2);
				
				% "Specifically, this class implements the pairwise or apparent transfer entropy, i.e. we compute the
				% transfer that appears to come from a single source variable, without examining any other potential sources"
				teCalc.setProperty('k_HISTORY', num2str(tau_history));
				teCalc.setProperty('k_TAU', num2str(tau_steps));
				teCalc.initialise(1, micro_dim, macro_dim);
				
				teCalc.setObservations(octaveToJavaDoubleMatrix(micro), ...
					octaveToJavaDoubleMatrix(macro));
				
				DD_nats = teCalc.computeAverageLocalOfObservations();
				DD = DD_nats/(1/log(2));
				
				all_DD.([fieldnames_micro{i} '_' fieldnames_macro{j}]) = DD;
				
			elseif strcmp(method, 'Discrete')
				
				% Class TransferEntropyCalculatorDiscrete: http://lizier.me/joseph/software/jidt/javadocs/v1.5/
				% Joe Lizier: "For discrete, the key is that we convert any of the Y_t , X_t-1 or Y_t-1 which are multi-
				% variate into a univariate time series, by basically converting into a larger alphabet space, e.g.
				% 2x binary variables convert into a base-4 univariate. This is done simply by using the utility
				% computeCombinedValues method, and setting the base or alphabet appropriately in the constructor,
				% and then you just use the calculator as is."
				
				teCalc.initialise();
				
				% We need to construct the joint values of the dest and source before we pass them in,
				% and need to use the matrix conversion routine when calling from Matlab/Octave:
				mUtils = javaObject('infodynamics.utils.MatrixUtils');
				teCalc.addObservations(mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(micro), 12), ...
					mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(macro), 1));
				DD = teCalc.computeAverageLocalOfObservations()
				
				all_DD.([fieldnames_micro{i} '_' fieldnames_macro{j}]) = DD;
			end
		end 
	end 
end 
