Steps to preprocesing the data:

Step1: Run the 0_mouse_indi.sh, which will run mouse_indi.sge which will run mouse_indi.R, to separate the samples from each batch
Step2: Run 1_mouse_asym_acr_sample_merge.sge will run mouse_asym_acr_sample_merge.R, integrating all samples together
Step3: Run 2_mouse_asym_acr_sample_harmony.sge will run mouse_asym_acr_sample_harmony.R, do harmony processing

other scripts' meaning:
sample_pisition.R -> plot the position of each sample of each batch, to decide the coordinates to crop the sample

mouse_asym_acr_sample_rpca.R -> runing integration with rpcA, NOT recommanded, taking ages
mouse_asym_acr_batch.R -> running integrationg across batches not samples, obselete

Xenium_example.R -> example script from Seurat
mouse_asym_within_batch.R -> run one example on batch 695
compare_before_after_harmony.R -> compare if there is any difference before and after harmony
compare_before_after_harmony_batch95.R -> compare if there is any difference before and after harmony within batch 695
