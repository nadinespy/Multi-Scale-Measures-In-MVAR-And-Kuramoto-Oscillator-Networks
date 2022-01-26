% INPUT:	 matrices for causal emergence for all parameters (each cell in will be a matrix of size A (number of different values for parameter 
%			#1), and B (number of different values for parameter #2)), using mmi & and micro #1-4; same for downward causation & causal 
%			decoupling;
%
% 			matrices for causal emergence for all parameters (each cell in will be a matrix of size A (number of different values for parameter 
%			#1), and B (number of different values for parameter #2)), using ccs & and micro #1-4; same for downward causation & causal 
%			decoupling; 
%
% 			matrices for practical causal emergence for all parameters (each cell in will be a matrix of size A (number of different values for 
%			parameter #1), and B (number of different values for parameter #2)), using macro #1 & and micro #1-4; 
%			same for downward causation & causal decoupling;  
%
% 			matrices for practical causal emergence for all parameters (each cell in will be a matrix of size A (number of different values for 
%			parameter #1), and B (number of different values for parameter #2)), using macro #2 & and micro #1-4; 
%			same for downward causation & causal decoupling. 
%
% OUTPUT:  ce_pract: struct with 16 fields denoting causal emergence, downward causation, & causal decoupling for all combinations of 
%			two macro and four micro variables.
%
%			ce_phiid_ccs & ce_phiid_mmi: structs with 12 fields each (one for mmi, one for ccs) denoting causal emergence, downward causation 
%			for all micro variables mmi & ccs, respectively.

% TODO: Find a way to input a flexible number of micro & macro variables

function [ce_pract] = store_pract_ce_in_struct_2x4(...
			ce_pract_macro1_micro1, ce_pract_macro1_micro2, ce_pract_macro1_micro3, ce_pract_macro1_micro4,...
			dc_pract_macro1_micro1, dc_pract_macro1_micro2, dc_pract_macro1_micro3, dc_pract_macro1_micro4, ...
			cd_pract_macro1_micro1, cd_pract_macro1_micro2, cd_pract_macro1_micro3, cd_pract_macro1_micro4,...
			ce_pract_macro2_micro1, ce_pract_macro2_micro2, ce_pract_macro2_micro3, ce_pract_macro2_micro4,...
			dc_pract_macro2_micro1, dc_pract_macro2_micro2, dc_pract_macro2_micro3, dc_pract_macro2_micro4,...
			cd_pract_macro2_micro1, cd_pract_macro2_micro2, cd_pract_macro2_micro3, cd_pract_macro2_micro4)
		
	% allocate three structs with three values each (causal emergence, downward causation, causal decoupling) 
	% - one for causal emergence using mmi, one for causal emergence using ccs, one for practical causal emergence
	ce_pract = [];
	
	% practical causal emergence using macro #1 and micro #1-4
	ce_pract.ce_pract_macro1_micro1 = ce_pract_macro1_micro1;
	ce_pract.ce_pract_macro1_micro2 = ce_pract_macro1_micro2;
	ce_pract.ce_pract_macro1_micro3 = ce_pract_macro1_micro3;
	ce_pract.ce_pract_macro1_micro4 = ce_pract_macro1_micro4;

	% downward causation using macro #1 and micro #1-4
	ce_pract.dc_pract_macro1_micro1 = dc_pract_macro1_micro1;
	ce_pract.dc_pract_macro1_micro2 = dc_pract_macro1_micro2;
	ce_pract.dc_pract_macro1_micro3 = dc_pract_macro1_micro3;
	ce_pract.dc_pract_macro1_micro4 = dc_pract_macro1_micro4;
	
	% causal decoupling using macro #1 and micro #1-4
	ce_pract.cd_pract_macro1_micro1 = cd_pract_macro1_micro1;
	ce_pract.cd_pract_macro1_micro2 = cd_pract_macro1_micro2;
	ce_pract.cd_pract_macro1_micro3 = cd_pract_macro1_micro3;
	ce_pract.cd_pract_macro1_micro4 = cd_pract_macro1_micro4;
	
	% practical causal emergence using macro #2 and micro #1-4
	ce_pract.ce_pract_macro2_micro1 = ce_pract_macro2_micro1;
	ce_pract.ce_pract_macro2_micro2 = ce_pract_macro2_micro2;
	ce_pract.ce_pract_macro2_micro3 = ce_pract_macro2_micro3;
	ce_pract.ce_pract_macro2_micro4 = ce_pract_macro2_micro4;

	% downward causation using macro #2 and micro #1-4
	ce_pract.dc_pract_macro2_micro1 = dc_pract_macro2_micro1;
	ce_pract.dc_pract_macro2_micro2 = dc_pract_macro2_micro2;
	ce_pract.dc_pract_macro2_micro3 = dc_pract_macro2_micro3;
	ce_pract.dc_pract_macro2_micro4 = dc_pract_macro2_micro4;
	
	% causal decoupling using macro #2 and micro #1-4
	ce_pract.cd_pract_macro2_micro1 = cd_pract_macro2_micro1;
	ce_pract.cd_pract_macro2_micro2 = cd_pract_macro2_micro2;
	ce_pract.cd_pract_macro2_micro3 = cd_pract_macro2_micro3;
	ce_pract.cd_pract_macro2_micro4 = cd_pract_macro2_micro4;
	
	% store input names in variables 1-48 (inputname(n) returns the input variable name of the nth argument 
	% as a string
	variable1 = inputname(1);
	variable2 = inputname(2);
	variable3 = inputname(3);
	variable4 = inputname(4);
	variable5 = inputname(5);
	variable6 = inputname(6);
	variable7 = inputname(7);
	variable8 = inputname(8);
	variable9 = inputname(9);
	variable10 = inputname(10);
	variable11 = inputname(11);
	variable12 = inputname(12);
	variable13 = inputname(13);
	variable14 = inputname(14);
	variable15 = inputname(15);
	variable16 = inputname(16);
	variable17 = inputname(17);
	variable18 = inputname(18);
	variable19 = inputname(19);
	variable20 = inputname(20);
	variable21 = inputname(21);
	variable22 = inputname(22);
	variable23 = inputname(23);
	variable24 = inputname(24);

	% change fieldnames to input names (by adding new fields and deleting the old ones)

	% practical causal emergence using macro #1 and micro #1-4
	ce_pract.(variable1) = ce_pract.ce_pract_macro1_micro1;
	ce_pract = rmfield(ce_pract,'ce_pract_macro1_micro1');
	ce_pract.(variable2) = ce_pract.ce_pract_macro1_micro2;
	ce_pract = rmfield(ce_pract,'ce_pract_macro1_micro2');
	ce_pract.(variable3) = ce_pract.ce_pract_macro1_micro3;
	ce_pract = rmfield(ce_pract,'ce_pract_macro1_micro3');
	ce_pract.(variable4) = ce_pract.ce_pract_macro1_micro4;
	ce_pract = rmfield(ce_pract,'ce_pract_macro1_micro4');

	% downward causation using macro #1 and micro #1-4
	ce_pract.(variable5) = ce_pract.dc_pract_macro1_micro1;
	ce_pract = rmfield(ce_pract,'dc_pract_macro1_micro1');
	ce_pract.(variable6) = ce_pract.dc_pract_macro1_micro2;
	ce_pract = rmfield(ce_pract,'dc_pract_macro1_micro2');
	ce_pract.(variable7) = ce_pract.dc_pract_macro1_micro3;
	ce_pract = rmfield(ce_pract,'dc_pract_macro1_micro3');
	ce_pract.(variable8) = ce_pract.dc_pract_macro1_micro4;
	ce_pract = rmfield(ce_pract,'dc_pract_macro1_micro4');
	
	% causal decoupling using macro #1 and micro #1-4
	ce_pract.(variable9) = ce_pract.cd_pract_macro1_micro1;
	ce_pract = rmfield(ce_pract,'cd_pract_macro1_micro1');
	ce_pract.(variable10) = ce_pract.cd_pract_macro1_micro2;
	ce_pract = rmfield(ce_pract,'cd_pract_macro1_micro2');
	ce_pract.(variable11) = ce_pract.cd_pract_macro1_micro3;
	ce_pract = rmfield(ce_pract,'cd_pract_macro1_micro3');
	ce_pract.(variable12) = ce_pract.cd_pract_macro1_micro4;
	ce_pract = rmfield(ce_pract,'cd_pract_macro1_micro4');
	
	% practical causal emergence using macro #2 and micro #1-4
	ce_pract.(variable13) = ce_pract.ce_pract_macro2_micro1;
	ce_pract = rmfield(ce_pract,'ce_pract_macro2_micro1');
	ce_pract.(variable14) = ce_pract.ce_pract_macro2_micro2;
	ce_pract = rmfield(ce_pract,'ce_pract_macro2_micro2');
	ce_pract.(variable15) = ce_pract.ce_pract_macro2_micro3;
	ce_pract = rmfield(ce_pract,'ce_pract_macro2_micro3');
	ce_pract.(variable16) = ce_pract.ce_pract_macro2_micro4;
	ce_pract = rmfield(ce_pract,'ce_pract_macro2_micro4');

	% downward causation using macro #2 and micro #1-4
	ce_pract.(variable17) = ce_pract.dc_pract_macro2_micro1;
	ce_pract = rmfield(ce_pract,'dc_pract_macro2_micro1');
	ce_pract.(variable18) = ce_pract.dc_pract_macro2_micro2;
	ce_pract = rmfield(ce_pract,'dc_pract_macro2_micro2');
	ce_pract.(variable19) = ce_pract.dc_pract_macro2_micro3;
	ce_pract = rmfield(ce_pract,'dc_pract_macro2_micro3');
	ce_pract.(variable20) = ce_pract.dc_pract_macro2_micro4;
	ce_pract = rmfield(ce_pract,'dc_pract_macro2_micro4');
	
	% causal decoupling using macro #2 and micro #1-4
	ce_pract.(variable21) = ce_pract.cd_pract_macro2_micro1;
	ce_pract = rmfield(ce_pract,'cd_pract_macro2_micro1');
	ce_pract.(variable22) = ce_pract.cd_pract_macro2_micro2;
	ce_pract = rmfield(ce_pract,'cd_pract_macro2_micro2');
	ce_pract.(variable23) = ce_pract.cd_pract_macro2_micro3;
	ce_pract = rmfield(ce_pract,'cd_pract_macro2_micro3');
	ce_pract.(variable24) = ce_pract.cd_pract_macro2_micro4;
	ce_pract = rmfield(ce_pract,'cd_pract_macro2_micro4');

end
