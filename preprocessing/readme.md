Steps to preprocesing the data:

Step1: Run the 0_mouse_indi.sh, which will run mouse_indi.sge which will run mouse_indi.R, to separate the samples from each batch
Step2: Run 1_mouse_asym_acr_sample_merge.sge will run mouse_asym_acr_sample_merge.R, integrating all samples together
Step3: Run 2_mouse_asym_acr_sample_harmony.sge will run mouse_asym_acr_sample_harmony.R, do harmony processing
