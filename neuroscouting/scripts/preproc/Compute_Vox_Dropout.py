from __future__ import division
import numpy as np
import utils as utils
import sys
import os


# Determine if the coverage of the functional scans results in needing to drop certain nodes due to sub-threshold coverage. 
#	a. Scripts: [1] /scripts/preproc/Compute_Vox_Dropout.py [2] /scripts/preproc/utils.py --> functions load_nifti3D, ROIcheck
#	b. Must edit contents of script for atlas-specific names, numbers, etc. 
#		i. [Compute_Vox_Dropout.py] Change 'nNodes', 'label', and remove the '0/1' used for debugging purposes
#
#==> If there are any missing nodes, the script with save a numpy file with...
#==>==> The index numbers for the missing nodes across all subjects. If x>0 subjects are missing 1 node, that node is considered missing for the entire group
#==>==> The index numbers for the common nodes. This is a list of all the nodes that had coverage in all subjects. 


### Assess which ROIs don't have enough voxels as per a set threshold.
### Requires an atlas/parcellation in a single volume (i.e. each node has a diff value) that is warped to the subject's T1W sapce and is the same resolution. 


##################################
### Atlases used in this study ###
################################## 

### Brainnetome ###

### Threshold Probability .25 (Cortex, Subcortex):
atlas_str = 'brainnetome/BNA_MPM_thr25_125mm_funcRes.nii.gz'

### DO NOT USE THIS ONE. IT HAS A 3-VOXEL NODE THAT DISAPPEARS AT FUNCTIONAL RESOLUTION 
#atlas_str = 'brainnetome/BN_Atlas_274_combined_funcRes.nii.gz'


### Schaefer 400 Node ###

### Yeo 7 Network (Cortex): 
#atlas_str = 'Schaefer2018_LocalGlobal/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm_funcRes.nii.gz'

### Yeo 17 Network (Cortex):
#atlas_str = 'Schaefer2018_LocalGlobal/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm_funcRes.nii.gz'

def main(argv = sys.argv):

	subDir = '/home/despoB/adameich/Experiments/nScout/data/derivatives/denoise_out'
	sublist = utils.get_sublist(subDir) 
	print(sublist)

	voxThresh = .25
	nNodes = 246 # 246 400
	label = 'Brainnetome' # 'Brainnetome' 'Schaefer400node'
	missingMat = np.zeros((len(sublist), 2)) * np.nan
	missingIdxsMat = []
	missingIdxsMatnSub = []
	unique_Missing_Nodes_raw = []

	for subIdx, sub in enumerate(sublist):
		missingIdxs_Per_Sub = []
		tmp = 0 
		for run in ['001', '002']:
			atlas = utils.load_nifti3D('/home/despoB/adameich/Experiments/nScout/data/resources/parcellations/%s' %(atlas_str))
			dataMask = utils.load_nifti3D('/home/despoB/adameich/Experiments/nScout/data/derivatives/fmriprep/sub-%s/func/sub-%s_task-rest_run-%s_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz' %(sub, sub, run))
			badrois,badcount,badroilist = utils.ROIcheck(atlas, dataMask, voxThresh)
			missingIdxsMat.append(badrois)
			missingIdxs_Per_Sub.append(badrois)
			unique_Missing_Nodes_raw = np.concatenate((unique_Missing_Nodes_raw, badroilist))
			
			missingMat[subIdx, tmp] = badcount / nNodes
			tmp+=1
				
	
	commonRoiMask = np.sum(missingIdxsMat,0)==0
	missingMatnRuns = np.sum(missingIdxsMat,0)

	unique_Missing_Nodes = np.unique(unique_Missing_Nodes_raw)
	unique_Missing_Nodes -= 1 #The badroilist from ROIcheck returns the LIST values as 1-400 to match the values in the atlas mask, but we want them in 0-399 for python. 

	commonROIs = [x for x in range(nNodes) if x not in  unique_Missing_Nodes]
	commonROIs = np.asarray(commonROIs) 
	
	#Save common ROI list
	np.save('/home/despoB/adameich/Experiments/nScout/data/analyses/commonROIs_list_thresh-%s_%s.npy' %(voxThresh, label), commonROIs)

	# Save missing ROI list
	np.save('/home/despoB/adameich/Experiments/nScout/data/analyses/missingROIs_list_thresh-%s_%s.npy' %(voxThresh,label), unique_Missing_Nodes)


### MAIN SCRIPT ###

if __name__ == '__main__':
    main()
