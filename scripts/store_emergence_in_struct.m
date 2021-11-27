function [emergence_ccs, emergence_mmi, emergence_practical] = store_emergence_in_struct(synergy_capacity_mmi, causal_decoupling_mmi, ...
		downward_causation_mmi, synergy_capacity_ccs, causal_decoupling_ccs, downward_causation_ccs, ...
		synergy_capacity_practical_linear, causal_decoupling_practical_linear, downward_causation_practical_linear)
	
	emergence_ccs = [];
	emergence_mmi = [];
	emergence_practical = [];
	
	emergence_ccs.synergy_capacity_ccs = synergy_capacity_ccs;
	emergence_ccs.causal_decoupling_ccs = causal_decoupling_ccs;
	emergence_ccs.downward_causation_ccs = downward_causation_ccs;
	
	emergence_mmi.synergy_capacity_mmi = synergy_capacity_mmi;
	emergence_mmi.causal_decoupling_mmi = causal_decoupling_mmi;
	emergence_mmi.downward_causation_mmi = downward_causation_mmi;
	
	emergence_practical.synergy_capacity_practical_linear = synergy_capacity_practical_linear;
	emergence_practical.causal_decoupling_practical_linear = causal_decoupling_practical_linear;
	emergence_practical.downward_causation_practical_linear = downward_causation_practical_linear;
	
end