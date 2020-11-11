for subject in 091 092 093 ; do #074 ... 
	python ROIExtractor_JRC.py \
	/home/despoB/adameich/Experiments/nScout/data/derivatives/denoise_out \
	-targetsuffix space-MNI152NLin2009cAsym_desc-residuals_variant-24p_Acc6_009_1_bold.nii.gz \
	-subject ${subject} \
	-roimask /home/despoB/adameich/Experiments/nScout/data/resources/parcellations/brainnetome/BNA_MPM_thr25_125mm_funcRes.nii.gz \
	-brainmask space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz \
	-atlasName brainnetome \
	-input_fmriprep /home/despoB/adameich/Experiments/nScout/data/derivatives/fmriprep \
	--nodatacheck

done
