function phiid_atoms = get_phiid_atoms(micro, method, time_lag, red_func)
	
	if strcmp(lower(method), 'gaussian');
		phiid_atoms = struct2array(PhiIDFull(micro, time_lag, red_func))';
	elseif strcmp(lower(method), 'discrete');
		phiid_atoms = struct2array(PhiIDFullDiscrete(micro, time_lag, red_func))';
	end 
	
end 