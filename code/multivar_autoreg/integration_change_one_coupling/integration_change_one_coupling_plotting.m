
BIG_array = [];
for i = 1:length(coupling_matrices);
	
	c_temp = c(i);
	
	for j = 1:length(all_rmi);
		
		rmi = all_rmi(j);
		temp_matrix = cell2mat(table2array(integration(i,j)));
		
		for k = 1:n_samples_noise_corrs;
		
			big_array = [c_temp, rmi, temp_matrix(:,k)];
			BIG_array = [BIG_array; big_array];
		
		end

	end 
	
	disp(i)
	
end 

figure;
scatter3(BIG_array(:,1), BIG_array(:,2), BIG_array(:,3), 2, ".", 'blue');
title({'integration across RMI and changes for one coupling: for rmi,', 'noise correlations were randomly sampled'});
xlabel('coupling of one node'); 
ylabel('RMI'); 
zlabel('total correlation');

location = string(strcat(pathout_plots, 'integration', '.png'));
exportgraphics(gcf, location);



