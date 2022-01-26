% INPUT: 	two structs with all atoms for all parameters (each cell in the struct will be a matrix of size A (number of different values for parameter 			#1), and B (number of different values for parameter #2)), the first referring to mmi as a redundancy function, the second referring to ccs.

% OUTPUT: 	matrices of size A (number of different values for parameter #1), and B (number of different values for parameter #2)), for 
% 		causal emergence, downward causation, & causal decoupling, and two different redundancy functions (first mmi, then ccs), respectively. 

function [ce_mmi, dc_mmi, cd_mmi, ce_ccs, ...
		dc_ccs, cd_ccs] = phiid_ce(all_atoms_mmi, all_atoms_ccs)
	
	ce_ccs = all_atoms_ccs.str + ...
		all_atoms_ccs.stx + all_atoms_ccs.sty + all_atoms_ccs.sts;
	
	dc_ccs = all_atoms_ccs.str + all_atoms_ccs.stx + all_atoms_ccs.sty;
	
	ce_mmi = all_atoms_mmi.str + ...
		all_atoms_mmi.stx + all_atoms_mmi.sty + all_atoms_mmi.sts;
	
	dc_mmi = all_atoms_mmi.str + all_atoms_mmi.stx + all_atoms_mmi.sty;
	
	cd_ccs = ce_ccs - dc_ccs;
	cd_mmi = ce_mmi - dc_mmi;
	
end 
