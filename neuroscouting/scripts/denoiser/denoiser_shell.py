import os, sys
import numpy as np

runs = ['001', '002']
epi_space = 'MNI152NLin2009cAsym' # 'T1w'

for sub in [93]: #np.arange(1,73):
	if sub in [8, 19, 56]:
		continue
	for run in runs:

		img_file_str = '/home/despoB/adameich/Experiments/nScout/data/derivatives/fmriprep/sub-%03d/func/sub-%03d_task-rest_run-%s_space-%s_desc-preproc_bold.nii.gz' %(sub, sub, run, epi_space)
		tsv_file_str = '/home/despoB/adameich/Experiments/nScout/data/derivatives/fmriprep/sub-%03d/func/sub-%03d_task-rest_run-%s_desc-confounds_regressors.tsv' %(sub, sub, run)
		out_path = '/home/despoB/adameich/Experiments/nScout/data/derivatives/denoise_out/sub-%03d' %(sub)

		if not os.path.isdir(out_path):
			os.mkdir(out_path)

		hp_filter = 0.009
		lp_filter = 0.1
		BASE_CONFOUNDS ='a_comp_cor_00 a_comp_cor_01 a_comp_cor_02 a_comp_cor_03 a_comp_cor_04 a_comp_cor_05 trans_x trans_x_derivative1 trans_x_derivative1_power2 trans_x_power2 trans_y trans_y_derivative1 trans_y_power2 trans_y_derivative1_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_power2 rot_x_derivative1_power2 rot_y rot_y_derivative1 rot_y_power2 rot_y_derivative1_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2'

		GSR_COLS = 'global_signal global_signal_derivative1 global_signal_derivative1_power2 global_signal_power2'

		col_names_str = BASE_CONFOUNDS # + ' ' + GSR_COLS
		fd_col_name_str = 'framewise_displacement'
		FD_thr = 0.3
		bids = True
		template_str = '/home/despoB/adameich/Experiments/nScout/scripts/denoiser/report_template.html'
		strategy_name_str = '24p_Acc6_009_1' # '24p_Acc6_009_1' '24p_Acc6_009_1_GSR4'

		system_call = 'python run_denoiser_test.py %s %s %s --hp_filter %s --lp_filter %s --col_names %s --fd_col_name %s --FD_thr %s --bids %s --strategy_name %s --template_file %s' %(img_file_str, tsv_file_str, out_path, hp_filter, lp_filter, col_names_str, fd_col_name_str, FD_thr, bids, strategy_name_str, template_str)
		print(system_call)



		os.system(system_call)
