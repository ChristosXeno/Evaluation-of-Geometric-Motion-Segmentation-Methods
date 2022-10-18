# Evaluation-of-Geometric-Motion-Segmentation-Methods
This repository contains code that was used to gather and the results for my 'Evaluation of Geometric Motion Segmentation Methods' thesis

In order to run a segmentation method on a particular dataset open the 'Thesis Motion Segmentation Methods' folder and then select the folder of the method of choice, from there, select the main script that is to be run for this and then alter the 'datasets_folder' and 'dataset_folder' variable so that the appropriate file path to the datasets folder is selected as well as the wanted dataset. From there, run the program from this main file. A list of the folders, the segmentation algorithms in which they run as well as their main script are listed below.
Format: Method = Folder Name (main script)
- RSIM = Robust-shape-interaction-matrix-master (main.m)
- LRSSC = greedysc_codes (run_sc_motionseg.m)
- GPCA, LSA and RANSAC = HopkinsMultiviewMultibodyCode (MoSegTest.m)
- T-Linkage = TLinkage (demo.m)
- RPA = RPA_v02 (demo.m)

To add a method into the program simply import it to the 'Thesis Motion Segmentation Methods' folder and ensure that the 'alter_data.m' script is ran after the appropriate dataset is imported for use in this code.

A new dataset can be added to this script as long as it contains information pertaining to the table in figure 3.1. If this exists, then it can be added and a new condition must then be created in the 'alter_data.m' script in order to allow it to be used within the program.

Note: The segmentation methods do not need to be ran in order for their results to be found as they are already saved for all scenarios described in the thesis report. The 'ResultsProcessing.m' file can be used to process and plot these.
