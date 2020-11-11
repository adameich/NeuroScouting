from __future__ import division
import numpy as np
import sys
import os
from nilearn.image import resample_to_img


### Resamples an atlas (i.e. source image) in MNI space into the resolution of the functional scans collected (3.0 x 3.0 x 3.0). 
### This does not perform any registration, simply a resampling, so it doesn't matter WHICH epi is used.

##################################
### Atlases used in this study ###
################################## 

### Brainnetome ###

### Threshold Probability .25 (Cortex, Subcortex):
#atlas_str = 'brainnetome/BNA_MPM_thr25_125mm.nii.gz'

### ([I assume] Threshold Probability .25) with Cerebellum (Cortex, Subcortex, Cerebellum):
#atlas_str = 'brainnetome/BN_Atlas_274_combined.nii.gz'

### Schaefer 400 Node ###

### Yeo 7 Network (Cortex): 
#atlas_str = 'Schaefer2018_LocalGlobal/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'

### Yeo 17 Network (Cortex):
atlas_str = 'Schaefer2018_LocalGlobal/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii.gz'



sourceImage = '/home/despoB/adameich/Experiments/nScout/data/resources/parcellations/%s' %(atlas_str)

refImage = '/home/despoB/adameich/Experiments/nScout/data/derivatives/fmriprep/sub-001/func/sub-001_task-rest_run-001_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'

print('\nResampling %s atlas [%s] in MNI space into epi resolution [3.0 x 3.0 x 3.0]' %(atlas_str.split('/')[0], atlas_str.split('/')[1]))

resampledMask = resample_to_img(sourceImage, refImage, interpolation='nearest')

resampledMask.to_filename('/home/despoB/adameich/Experiments/nScout/data/resources/parcellations/%s/%s_funcRes.nii.gz' %(atlas_str.split('/')[0], atlas_str.split('/')[1].split('.')[0]))
