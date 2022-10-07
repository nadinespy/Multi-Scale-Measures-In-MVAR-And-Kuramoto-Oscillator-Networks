function all_DD = get_DD(micro_variables, macro_variables, method, time_lag, time_step, varargin)

	% Function description: get_DD() calculates DD for any micro and 
	% macro variable specified in [micro_variables] and [macro_variables],
	% and the parameters given; 
	
	% Inputs:	
	%
	% Required:	micro_variables			1x1 struct with micro 
	%							variable names in fields, 
	%							and micro variable 
	%							time-series as values
	%		micro_variables			1x1 struct with macro 
	%							variable names in fields, 
	%							and macro variable 
	%							time-series as values
	%		method				character array
	%		time_lag				double
	%		time_step				double
	%
	% Optional: kraskov_param			double
	%
	% Outputs:	all_DD				1x1 struct with with micro- 
	%							macro combinations in  
	%							fields, and DD as values
	
	% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'micro_variables', @isstruct);
	addRequired(p,'macro_variables', @isstruct);
	addRequired(p,'method', @ischar);
	addRequired(p,'time_lag', @isdouble);
	addRequired(p,'time_step', @isdouble);
	
	% optional name-value pair variables: 
	default_kraskov_param = 3;
	addParameter(p,'kraskov_param', default_kraskov_param, @isdouble);
	
	parse(p, micro_variables, macro_variables, method, time_lag, ...
		time_step, varargin{:});

	time_lag				= p.Results.time_lag;
	time_step				= p.Results.time_step;
	micro_variables			= p.Results.micro_variables;
	macro_variables			= p.Results.macro_variables;
	kraskov_param			= p.Results.kraskov_param;

	% get number of micro and macro variables, respectively
	n_micro_variables = length(fieldnames(micro_variables)); 
	n_macro_variables = length(fieldnames(macro_variables)); 
	
	% get names of micro and macro variables, respectively
	fieldnames_macro = fieldnames(macro_variables);
	fieldnames_micro = fieldnames(micro_variables);
	
	if strcmp(method, 'Kraskov')
		teCalc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorMultiVariateKraskov');
		teCalc.setProperty('k', num2str(kraskov_param));
	elseif strcmp(method, 'Gaussian');
		teCalc = javaObject('infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorMultiVariateGaussian');
	end 

	% loop over all micro and macro variables and calculate practical CE
	for i = 1:n_micro_variables;
		micro = micro_variables.(fieldnames_micro{i});
		
		for j = 1:n_macro_variables;
			macro = macro_variables.(fieldnames_macro{j});
			
			micro_dim = size(micro, 1);
			macro_dim = size(macro, 1);
	
			if strcmp(lower(method), 'kraskov') | strcmp(lower(method), 'gaussian');

				
				% "Specifically, this class implements the pairwise
				% or apparent transfer entropy, i.e. we compute the
				% transfer that appears to come from a single source 
				% variable, without examining any other potential sources"
				teCalc.setProperty('k_HISTORY', num2str(time_lag));
				teCalc.setProperty('k_TAU', num2str(time_step));
				teCalc.initialise(1, micro_dim, macro_dim);
				
				teCalc.setObservations(octaveToJavaDoubleMatrix(micro'), ...
					octaveToJavaDoubleMatrix(macro'));
				
				DD_nats = teCalc.computeAverageLocalOfObservations();
				DD = DD_nats/(1/log(2));
				
				all_DD.([fieldnames_micro{i} '_' fieldnames_macro{j}]) = DD;
				
			elseif strcmp(lower(method), 'discrete');
				
				% we need the alphabet size / number of states for 
				% each sample of the source (micro) and target (macro); 
				% for a multivariate variable, this means the number 
				% of states for each joint variable
				
				n_micro_states = numel(unique(micro));			% number of micro states
				n_macro_states = numel(unique(macro));			% number of macro states
				n_joint_micro_states = n_micro_states^size(micro,1);	% number ofjoint micro states
				
				% Joe Lizier on the problem of a too large state space: 
				% "The state space of joining[, e. g.,] 256 binary 
				% variables is just too large: 2^256 = 1.157920892×10⁷⁷. 
				% The estimator will fall over in trying to allocate 
				% memory to count each possible joint sample here, and 
				% whilst I do have code coming that will run the 
				% estimation without allocating such space it still 
				% won't work properly because that space is way too 
				% large for you to ever have enough samples to estimate 
				% properly. Roughly speaking, your number of samples 
				% should be 3x (minimum) or 10x (better) the number of 
				% joint states that you're likely to see.
				
				% TE estimator is built using number of joint micro states
				
				% class TransferEntropyCalculatorDiscrete: 
				% http://lizier.me/joseph/software/jidt/javadocs/v1.5/
				%
				% Joe Lizier: "For discrete, the key is that we convert 
				% any of the Y_t , X_t-1 or Y_t-1 which are multi-
				% variate into a univariate time series, by basically 
				% converting into a larger alphabet space, e.g.
				% 2x binary variables convert into a base-4 univariate. 
				% This is done simply by using the utility
				% computeCombinedValues method, and setting the base 
				% or alphabet appropriately in the constructor,
				% and then you just use the calculator as is."
				
				teCalc = javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete', ...
					n_joint_micro_states, time_lag);

				teCalc.initialise();

				% we need to construct the joint values of the dest 
				% and source before we pass them in, so 
				% mUtils.computeCombinedValues needs the number of 
				% states of micro/macro, respectively
				mUtils = javaObject('infodynamics.utils.MatrixUtils');
				
				% no need to use mUtils.computeCombinedValues(), if macro variable has only dimensions
				if size(macro, 1) > 1
					teCalc.addObservations(mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(micro'), ...
					n_micro_states), mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(macro'), ...
					n_macro_states));
				else 
					teCalc.addObservations(mUtils.computeCombinedValues(octaveToJavaDoubleMatrix(micro'), ...
					n_micro_states), macro');
				end 
				
				DD = teCalc.computeAverageLocalOfObservations();
				
				all_DD.([fieldnames_micro{i} '_' fieldnames_macro{j}]) = DD;
			end
		end 
	end 
end 
