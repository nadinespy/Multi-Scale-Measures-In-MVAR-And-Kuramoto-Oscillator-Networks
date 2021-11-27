function [all_atoms_mmi, all_atoms_ccs] = store_atoms_in_struct(phiid_mmi, phiid_ccs);
	all_atoms_mmi = [];
	
	all_atoms_mmi.rtr = squeeze(phiid_mmi(1,:,:));
	all_atoms_mmi.rtx = squeeze(phiid_mmi(2,:,:));
	all_atoms_mmi.rty = squeeze(phiid_mmi(3,:,:));
	all_atoms_mmi.rts = squeeze(phiid_mmi(4,:,:));
	all_atoms_mmi.xtr = squeeze(phiid_mmi(5,:,:));
	all_atoms_mmi.xtx = squeeze(phiid_mmi(6,:,:));
	all_atoms_mmi.xty = squeeze(phiid_mmi(7,:,:));
	all_atoms_mmi.xts = squeeze(phiid_mmi(8,:,:));
	all_atoms_mmi.ytr = squeeze(phiid_mmi(9,:,:));
	all_atoms_mmi.ytx = squeeze(phiid_mmi(10,:,:));
	all_atoms_mmi.yty = squeeze(phiid_mmi(11,:,:));
	all_atoms_mmi.yts = squeeze(phiid_mmi(12,:,:));
	all_atoms_mmi.str = squeeze(phiid_mmi(13,:,:));
	all_atoms_mmi.stx = squeeze(phiid_mmi(14,:,:));
	all_atoms_mmi.sty = squeeze(phiid_mmi(15,:,:));
	all_atoms_mmi.sts = squeeze(phiid_mmi(16,:,:));
	
	all_atoms_ccs = [];
	all_atoms_ccs.rtr = squeeze(phiid_ccs(1,:,:));
	all_atoms_ccs.rtx = squeeze(phiid_ccs(2,:,:));
	all_atoms_ccs.rty = squeeze(phiid_ccs(3,:,:));
	all_atoms_ccs.rts = squeeze(phiid_ccs(4,:,:));
	all_atoms_ccs.xtr = squeeze(phiid_ccs(5,:,:));
	all_atoms_ccs.xtx = squeeze(phiid_ccs(6,:,:));
	all_atoms_ccs.xty = squeeze(phiid_ccs(7,:,:));
	all_atoms_ccs.xts = squeeze(phiid_ccs(8,:,:));
	all_atoms_ccs.ytr = squeeze(phiid_ccs(9,:,:));
	all_atoms_ccs.ytx = squeeze(phiid_ccs(10,:,:));
	all_atoms_ccs.yty = squeeze(phiid_ccs(11,:,:));
	all_atoms_ccs.yts = squeeze(phiid_ccs(12,:,:));
	all_atoms_ccs.str = squeeze(phiid_ccs(13,:,:));
	all_atoms_ccs.stx = squeeze(phiid_ccs(14,:,:));
	all_atoms_ccs.sty = squeeze(phiid_ccs(15,:,:));
	all_atoms_ccs.sts = squeeze(phiid_ccs(16,:,:));
	
end 