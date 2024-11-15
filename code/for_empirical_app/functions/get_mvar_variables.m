function get_mvar_variables(network, model_sim_params, measure_params, ...
		coupling_matrices, pathout);
% get_mvar_variables() - generates micro and macro variables in MVAR networks.
% 
% Takes as inputs, amongst others, arrays with measure parameter values common to 
% all emergence calculations ([measure_params], required), model parameters used 
% for simulation ([model_sim_params], required), and coupling matrices 
% (coupling_matrices, required). It then loops over time-lengths and model parameters.
%
% Example: get_mvar_variables(network, model_sim_params, measure_params, ...
%          coupling_matrices, pathout)
% 
% INPUTS - required:	
%    network                -             char array
%    model_sim_params	    -		      1x1 struct with fields 
%							A (vector with doubles), beta, 
%							(vector with doubles), 	
%							intra_comm_size (int),
%							n_communities (int)
%
%    measure_params         -             1x1 struct with fields
%							'measures', 'methods',
%							'time_lags', 'time_lengths',
%							'kraskov_params', 'disc_methods',
%							' bins'
%
%							'measures': cell array with chars
%							'methods': cell array with chars
%							'time_lags': double array
%							'time_lengths': double array
%							'kraskov_params': double array
%							'disc_methods': cell arrays with chars
%							'bins': int array
%
%    coupling_matrices      -             3D double matrix of size 
%							[intra_comm_size*n_communities
%							 x intra_comm_size*n_communities
%							 x length(A)]
%
%    pathout                -             1x2 struct with field 
%							indicating path to output 
%							for simulated time-series, 
%							and synchronies
%
% OUTPUTS: 
%    micro and macro        -			saved in pathout
%    variables	 
%

% use inputParser to declare required & optional variables
	p = inputParser;
	
	% required variables
	addRequired(p,'network', @ischar);
	addRequired(p,'model_sim_params', @isstruct);
	addRequired(p,'measure_params', @isstruct);
	addRequired(p,'coupling_matrices', @isdouble);
	addRequired(p,'pathout', @isstruct);
	
	parse(p, network, model_sim_params, measure_params, coupling_matrices, pathout);
	
	network				= p.Results.network;
	model_sim_params			= p.Results.model_sim_params;
	measure_params			= p.Results.measure_params;
	coupling_matrices			= p.Results.coupling_matrices;
	pathout				= p.Results.pathout;
	
	% extract cell arrays from structs
	model_sim_params_fieldnames	= fieldnames(model_sim_params);
	pathout_fieldnames		= fieldnames(pathout);
	
	pathout				= pathout.(pathout_fieldnames{1});
	model_params1			= model_sim_params.(model_sim_params_fieldnames{1});
	model_params2			= model_sim_params.(model_sim_params_fieldnames{2});
	time_lengths			= measure_params.time_lengths;
	time_lag				= model_sim_params.time_lag_for_model;

	% SIMULATE MVAR NETWORK: outputs of sim_mvar_network() are
	%	- time-series of nodes				(MICRO)

	% GET FURTHER MICRO AND MACRO VARIABLES:
	%	- components of phases (half # of nodes,  (MICRO)
	%       if # of nodes < 2, otherwise # of nodes)  
	%	- summation over nodes				(MACRO)
	%	- summation over the exponent of nodes	(MACRO)

	for q = 1:length(time_lengths);
		time_length_str = num2str(time_lengths(q));
	
		for i = 1:size(coupling_matrices, 3);
			coupling_matrix = coupling_matrices(:,:,i);
			model_param1_str = param2str(model_params1(i));
		
			for j = 1:length(model_params2)
				model_param2_str = param2str(model_params2(j));
				
				fprintf('get_mvar_variables - loop indices: time_series_length: %d, model_param1: %d, model_param2: %d\n', ...
					q, i, j);
				
				rng(1);
				% simulate km oscillators using Shanahan's code
				nodes = sim_mvar_network(time_lengths(q), model_params2(j), coupling_matrix, time_lag);
				
				% MICRO VARIABLES:
				%	- time-series of nodes
				%	- components of phases
				
				% do reconstruction ICA with nodes: first get weights for each variable and each feature 
				% (output of rica() will be a matrix of size [size(nodes,2) * number of features], then 
				% multiply this matrix with input matrix to get time-series of independent 
				% components/projection of each data point in the component space
				
				if size(nodes, 1) <= 2
					n_features = size(nodes,1);
				else
					n_features = size(nodes,1)/2;
				end
				
				reconstruction_ica = rica(nodes', n_features,'IterationLimit',100);
				rica_nodes = (nodes' * reconstruction_ica.TransformWeights)';
				S.(['rica' num2str(n_features) '_nodes']) = rica_nodes;
				
				% MACRO VARIABLES:
				%	- summation over nodes
				%	- summation over the exponent of nodes
				
				% get sum of nodes
				sum_nodes = zeros(1, time_lengths(q));
				for k = 1:(size(nodes,1));
					sum_nodes = sum_nodes + nodes(k,:);
				end
				
				% get sum of exponent of nodes
				sum_exp_nodes = zeros(1, time_lengths(q));
				for k = 1:(size(nodes,1));
					sum_exp_nodes = sum_exp_nodes + exp(nodes(k,:));
				end
				
				% variable names consist of  network name + variable name + value of coupling + 
				% value of corr err + number of datapoints + time lag
				save([pathout network '_nodes_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'nodes');
				save([pathout network '_rica' num2str(n_features) '_nodes_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'-struct', 'S');
				save([pathout network '_sum_nodes_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'sum_nodes');
				save([pathout network '_sum_exp_nodes_' model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
					'sum_exp_nodes');
				
			end

			clear S;
		end
	end
end