% downsampling particular variables

ns = 5000;
dsf = 2;

for q = 1:length(time_lengths);
	time_length_str = num2str(time_lengths(q));
						   
	for i = 1:size(coupling_matrices, 3);
		coupling_matrix = coupling_matrices(:,:,i);
		model_param1_str = param2str(A(i));
							   
		for j = 1:length(beta)
			model_param2_str = param2str(beta(j));
						
			macro_variable = struct2array(load([pathout_data_sim_time_series ...
				network '_' 'full_system_sync' '_' model_param1_str '_' ...
				model_param2_str '_' time_length_str '.mat']));
								   
			m = floor((length(macro_variable))/dsf);
			X = zeros(size(macro_variable, 2),m);
			X = downsample(macro_variable,dsf);
								   
			full_system_sync = X;				   
			save([pathout_data_sim_time_series network '_full_system_sync_' ...
				model_param1_str '_' model_param2_str '_' time_length_str '.mat'], ...
				'full_system_sync');
			
			clear X
			clear full_system_sync
		end
	end
end
					   