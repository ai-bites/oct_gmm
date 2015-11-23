This is an outlier detection approach to OCT data classification using Gaussian Mixture Model. 

To set it up to run for your data, follow the below steps:

1. Get the code from the repo.
2. Create folders, 'data' and 'duke_data' that contain the dataset volumes in *.mat file format. These folders in turn contain folders 'DME' and 'normal' under them, which contain the corresponding volumes.
3. Add the folders 'data', 'duke_data' and 'utils' to Matlab path.
4. If you are running the project for the first time, set `do_preprocess` to 1, which runs NL Means on each frame in each volume and stores the result. To just use the saved result this can then be set to 0.
5. You are now all set to run the files for different cases:
	* For SERI dataset + intensity features run: run.m
	* For Duke dataset + intensity features run: run_duke.m
	* For SERI dataset + LBP features run: run_lbp.m
	* For Duke dataset + LBP features, run: run_lbp_duke.m
