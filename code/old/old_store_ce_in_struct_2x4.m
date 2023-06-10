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

function [ce_phiid_ccs, ce_phiid_mmi, ce_pract] = store_ce_in_struct_2x4(...
			ce_phiid_mmi_micro1, ce_phiid_mmi_micro2, ce_phiid_mmi_micro3, ce_phiid_mmi_micro4,...
			dc_phiid_mmi_micro1, dc_phiid_mmi_micro2, dc_phiid_mmi_micro3, dc_phiid_mmi_micro4, ...
			cd_phiid_mmi_micro1, cd_phiid_mmi_micro2, cd_phiid_mmi_micro3, cd_phiid_mmi_micro4,...
			ce_phiid_ccs_micro1, ce_phiid_ccs_micro2, ce_phiid_ccs_micro3, ce_phiid_ccs_micro4,...
			dc_phiid_ccs_micro1, dc_phiid_ccs_micro2, dc_phiid_ccs_micro3, dc_phiid_ccs_micro4, ...
			cd_phiid_ccs_micro1, cd_phiid_ccs_micro2, cd_phiid_ccs_micro3, cd_phiid_ccs_micro4,...
			ce_pract_macro1_micro1, ce_pract_macro1_micro2, ce_pract_macro1_micro3, ce_pract_macro1_micro4,...
			dc_pract_macro1_micro1, dc_pract_macro1_micro2, dc_pract_macro1_micro3, dc_pract_macro1_micro4, ...
			cd_pract_macro1_micro1, cd_pract_macro1_micro2, cd_pract_macro1_micro3, cd_pract_macro1_micro4,...
			ce_pract_macro2_micro1, ce_pract_macro2_micro2, ce_pract_macro2_micro3, ce_pract_macro2_micro4,...
			dc_pract_macro2_micro1, dc_pract_macro2_micro2, dc_pract_macro2_micro3, dc_pract_macro2_micro4,...
			cd_pract_macro2_micro1, cd_pract_macro2_micro2, cd_pract_macro2_micro3, cd_pract_macro2_micro4)
		
	% allocate three structs with three values each (causal emergence, downward causation, causal decoupling) 
	% - one for causal emergence using mmi, one for causal emergence using ccs, one for practical causal emergence
	ce_phiid_ccs = [];
	ce_phiid_mmi = [];
	ce_pract = [];
	
	% causal emergence for mmi using micro #1-4
	ce_phiid_mmi.ce_phiid_mmi_micro1 = ce_phiid_mmi_micro1;
	ce_phiid_mmi.ce_phiid_mmi_micro2 = ce_phiid_mmi_micro2;
	ce_phiid_mmi.ce_phiid_mmi_micro3 = ce_phiid_mmi_micro3;
	ce_phiid_mmi.ce_phiid_mmi_micro4 = ce_phiid_mmi_micro4;
	
	% downward causation for mmi using micro #1-4
	ce_phiid_mmi.dc_phiid_mmi_micro1 = dc_phiid_mmi_micro1;
	ce_phiid_mmi.dc_phiid_mmi_micro2 = dc_phiid_mmi_micro2;
	ce_phiid_mmi.dc_phiid_mmi_micro3 = dc_phiid_mmi_micro3;
	ce_phiid_mmi.dc_phiid_mmi_micro4 = dc_phiid_mmi_micro4;
	
	% causal decoupling for mmi using micro #1-4
	ce_phiid_mmi.cd_phiid_mmi_micro1 = cd_phiid_mmi_micro1;
	ce_phiid_mmi.cd_phiid_mmi_micro2 = cd_phiid_mmi_micro2;
	ce_phiid_mmi.cd_phiid_mmi_micro3 = cd_phiid_mmi_micro3;
	ce_phiid_mmi.cd_phiid_mmi_micro4 = cd_phiid_mmi_micro4;
	
	% causal emergence for ccs using micro #1-4
	ce_phiid_ccs.ce_phiid_ccs_micro1 = ce_phiid_ccs_micro1;
	ce_phiid_ccs.ce_phiid_ccs_micro2 = ce_phiid_ccs_micro2;
	ce_phiid_ccs.ce_phiid_ccs_micro3 = ce_phiid_ccs_micro3;
	ce_phiid_ccs.ce_phiid_ccs_micro4 = ce_phiid_ccs_micro4;
	
	% downward causation for ccs using micro #1-4
	ce_phiid_ccs.dc_phiid_ccs_micro1 = dc_phiid_ccs_micro1;
	ce_phiid_ccs.dc_phiid_ccs_micro2 = dc_phiid_ccs_micro2;
	ce_phiid_ccs.dc_phiid_ccs_micro3 = dc_phiid_ccs_micro3;
	ce_phiid_ccs.dc_phiid_ccs_micro4 = dc_phiid_ccs_micro4;
	
	% causal decoupling for ccs using micro #1-4
	ce_phiid_ccs.cd_phiid_ccs_micro1 = cd_phiid_ccs_micro1;
	ce_phiid_ccs.cd_phiid_ccs_micro2 = cd_phiid_ccs_micro2;
	ce_phiid_ccs.cd_phiid_ccs_micro3 = cd_phiid_ccs_micro3;
	ce_phiid_ccs.cd_phiid_ccs_micro4 = cd_phiid_ccs_micro4;

	
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
	variable25 = inputname(25);
	variable26 = inputname(26);
	variable27 = inputname(27);
	variable28 = inputname(28);
	variable29 = inputname(29);
	variable30 = inputname(30);
	variable31 = inputname(31);
	variable32 = inputname(32);
	variable33 = inputname(33);
	variable34 = inputname(34);
	variable35 = inputname(35);
	variable36 = inputname(36);
	variable37 = inputname(37);
	variable38 = inputname(38);
	variable39 = inputname(39);
	variable40 = inputname(40);
	variable41 = inputname(41);
	variable42 = inputname(42);
	variable43 = inputname(43);
	variable44 = inputname(44);
	variable45 = inputname(45);
	variable46 = inputname(46);
	variable47 = inputname(47);
	variable48 = inputname(48);

	% change fieldnames to input names (by adding new fields and deleting the old ones)
	
	% causal emergence for mmi using micro #1-4
	ce_phiid_mmi.(variable1) = ce_phiid_mmi.ce_phiid_mmi_micro1;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'ce_phiid_mmi_micro1');
	ce_phiid_mmi.(variable2) = ce_phiid_mmi.ce_phiid_mmi_micro2;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'ce_phiid_mmi_micro2');
	ce_phiid_mmi.(variable3) = ce_phiid_mmi.ce_phiid_mmi_micro3;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'ce_phiid_mmi_micro3');
	ce_phiid_mmi.(variable4) = ce_phiid_mmi.ce_phiid_mmi_micro4;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'ce_phiid_mmi_micro4');
	
	% downward causation for mmi using micro #1-4
	ce_phiid_mmi.(variable5) = ce_phiid_mmi.dc_phiid_mmi_micro1;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'dc_phiid_mmi_micro1');
	ce_phiid_mmi.(variable6) = ce_phiid_mmi.dc_phiid_mmi_micro2;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'dc_phiid_mmi_micro2');
	ce_phiid_mmi.(variable7) = ce_phiid_mmi.dc_phiid_mmi_micro3;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'dc_phiid_mmi_micro3');
	ce_phiid_mmi.(variable8) = ce_phiid_mmi.dc_phiid_mmi_micro4;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'dc_phiid_mmi_micro4');
	
	% causal decoupling for mmi using micro #1-4d_pract_macro2_micro4
	ce_phiid_mmi.(variable9) = ce_phiid_mmi.cd_phiid_mmi_micro1;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'cd_phiid_mmi_micro1');
	ce_phiid_mmi.(variable10) = ce_phiid_mmi.cd_phiid_mmi_micro2;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'cd_phiid_mmi_micro2');
	ce_phiid_mmi.(variable11) = ce_phiid_mmi.cd_phiid_mmi_micro3;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'cd_phiid_mmi_micro3');
	ce_phiid_mmi.(variable12) = ce_phiid_mmi.cd_phiid_mmi_micro4;
	ce_phiid_mmi = rmfield(ce_phiid_mmi,'cd_phiid_mmi_micro4');
	
	% causal emergence for ccs using micro #1-4
	ce_phiid_ccs.(variable13) = ce_phiid_ccs.ce_phiid_ccs_micro1;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'ce_phiid_ccs_micro1');
	ce_phiid_ccs.(variable14) = ce_phiid_ccs.ce_phiid_ccs_micro2;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'ce_phiid_ccs_micro2');
	ce_phiid_ccs.(variable15) = ce_phiid_ccs.ce_phiid_ccs_micro3;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'ce_phiid_ccs_micro3');
	ce_phiid_ccs.(variable16) = ce_phiid_ccs.ce_phiid_ccs_micro4;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'ce_phiid_ccs_micro4');
	
	% downward causation for ccs using micro #1-4
	ce_phiid_ccs.(variable17) = ce_phiid_ccs.dc_phiid_ccs_micro1;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'dc_phiid_ccs_micro1');
	ce_phiid_ccs.(variable18) = ce_phiid_ccs.dc_phiid_ccs_micro2;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'dc_phiid_ccs_micro2');
	ce_phiid_ccs.(variable19) = ce_phiid_ccs.dc_phiid_ccs_micro3;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'dc_phiid_ccs_micro3');
	ce_phiid_ccs.(variable20) = ce_phiid_ccs.dc_phiid_ccs_micro4;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'dc_phiid_ccs_micro4');
	
	% causal decoupling for ccs using micro #1-4
	ce_phiid_ccs.(variable21) = ce_phiid_ccs.cd_phiid_ccs_micro1;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'cd_phiid_ccs_micro1');
	ce_phiid_ccs.(variable22) = ce_phiid_ccs.cd_phiid_ccs_micro2;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'cd_phiid_ccs_micro2');
	ce_phiid_ccs.(variable23) = ce_phiid_ccs.cd_phiid_ccs_micro3;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'cd_phiid_ccs_micro3');
	ce_phiid_ccs.(variable24) = ce_phiid_ccs.cd_phiid_ccs_micro4;
	ce_phiid_ccs = rmfield(ce_phiid_ccs,'cd_phiid_ccs_micro4');

	% practical causal emergence using macro #1 and micro #1-4
	ce_pract.(variable25) = ce_pract.ce_pract_macro1_micro1;
	ce_pract = rmfield(ce_pract,'ce_pract_macro1_micro1');
	ce_pract.(variable26) = ce_pract.ce_pract_macro1_micro2;
	ce_pract = rmfield(ce_pract,'ce_pract_macro1_micro2');
	ce_pract.(variable27) = ce_pract.ce_pract_macro1_micro3;
	ce_pract = rmfield(ce_pract,'ce_pract_macro1_micro3');
	ce_pract.(variable28) = ce_pract.ce_pract_macro1_micro4;
	ce_pract = rmfield(ce_pract,'ce_pract_macro1_micro4');

	% downward causation using macro #1 and micro #1-4
	ce_pract.(variable29) = ce_pract.dc_pract_macro1_micro1;
	ce_pract = rmfield(ce_pract,'dc_pract_macro1_micro1');
	ce_pract.(variable30) = ce_pract.dc_pract_macro1_micro2;
	ce_pract = rmfield(ce_pract,'dc_pract_macro1_micro2');
	ce_pract.(variable31) = ce_pract.dc_pract_macro1_micro3;
	ce_pract = rmfield(ce_pract,'dc_pract_macro1_micro3');
	ce_pract.(variable32) = ce_pract.dc_pract_macro1_micro4;
	ce_pract = rmfield(ce_pract,'dc_pract_macro1_micro4');
	
	% causal decoupling using macro #1 and micro #1-4
	ce_pract.(variable33) = ce_pract.cd_pract_macro1_micro1;
	ce_pract = rmfield(ce_pract,'cd_pract_macro1_micro1');
	ce_pract.(variable34) = ce_pract.cd_pract_macro1_micro2;
	ce_pract = rmfield(ce_pract,'cd_pract_macro1_micro2');
	ce_pract.(variable35) = ce_pract.cd_pract_macro1_micro3;
	ce_pract = rmfield(ce_pract,'cd_pract_macro1_micro3');
	ce_pract.(variable36) = ce_pract.cd_pract_macro1_micro4;
	ce_pract = rmfield(ce_pract,'cd_pract_macro1_micro4');
	
	% practical causal emergence using macro #2 and micro #1-4
	ce_pract.(variable37) = ce_pract.ce_pract_macro2_micro1;
	ce_pract = rmfield(ce_pract,'ce_pract_macro2_micro1');
	ce_pract.(variable38) = ce_pract.ce_pract_macro2_micro2;
	ce_pract = rmfield(ce_pract,'ce_pract_macro2_micro2');
	ce_pract.(variable39) = ce_pract.ce_pract_macro2_micro3;
	ce_pract = rmfield(ce_pract,'ce_pract_macro2_micro3');
	ce_pract.(variable40) = ce_pract.ce_pract_macro2_micro4;
	ce_pract = rmfield(ce_pract,'ce_pract_macro2_micro4');

	% downward causation using macro #2 and micro #1-4
	ce_pract.(variable41) = ce_pract.dc_pract_macro2_micro1;
	ce_pract = rmfield(ce_pract,'dc_pract_macro2_micro1');
	ce_pract.(variable42) = ce_pract.dc_pract_macro2_micro2;
	ce_pract = rmfield(ce_pract,'dc_pract_macro2_micro2');
	ce_pract.(variable43) = ce_pract.dc_pract_macro2_micro3;
	ce_pract = rmfield(ce_pract,'dc_pract_macro2_micro3');
	ce_pract.(variable44) = ce_pract.dc_pract_macro2_micro4;
	ce_pract = rmfield(ce_pract,'dc_pract_macro2_micro4');
	
	% causal decoupling using macro #2 and micro #1-4
	ce_pract.(variable45) = ce_pract.cd_pract_macro2_micro1;
	ce_pract = rmfield(ce_pract,'cd_pract_macro2_micro1');
	ce_pract.(variable46) = ce_pract.cd_pract_macro2_micro2;
	ce_pract = rmfield(ce_pract,'cd_pract_macro2_micro2');
	ce_pract.(variable47) = ce_pract.cd_pract_macro2_micro3;
	ce_pract = rmfield(ce_pract,'cd_pract_macro2_micro3');
	ce_pract.(variable48) = ce_pract.cd_pract_macro2_micro4;
	ce_pract = rmfield(ce_pract,'cd_pract_macro2_micro4');

end
