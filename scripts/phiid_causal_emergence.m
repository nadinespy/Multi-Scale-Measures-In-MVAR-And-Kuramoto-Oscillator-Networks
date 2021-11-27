function [synergy_capacity_mmi, downward_causation_mmi, causal_decoupling_mmi, synergy_capacity_ccs, ...
		downward_causation_ccs, causal_decoupling_ccs] = phiid_causal_emergence(all_atoms_mmi, all_atoms_ccs)
	
	synergy_capacity_ccs = all_atoms_ccs.str + ...
		all_atoms_ccs.stx + all_atoms_ccs.sty + all_atoms_ccs.sts;
	
	downward_causation_ccs = all_atoms_ccs.str + all_atoms_ccs.stx + all_atoms_ccs.sty;
	
	synergy_capacity_mmi = all_atoms_mmi.str + ...
		all_atoms_mmi.stx + all_atoms_mmi.sty + all_atoms_mmi.sts;
	
	downward_causation_mmi = all_atoms_mmi.str + all_atoms_mmi.stx + all_atoms_mmi.sty;
	
	causal_decoupling_ccs = synergy_capacity_ccs - downward_causation_ccs;
	causal_decoupling_mmi = synergy_capacity_mmi - downward_causation_mmi;
	
end 