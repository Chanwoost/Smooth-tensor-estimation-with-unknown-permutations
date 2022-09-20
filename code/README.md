`functions_blk.R`: main functions used in the paper.

* `Borda2`: Estimate a signal tensor and permutation from a noisy and incomplete data tensor using Borda count estimation method.
* `Spectral`: Estimate a permuted signal tensor from a noisy data tensor using spectral method, which performs universal singualr value thresholding on the unfolded tensor.
* `LSE`: Estimate a permuted signal tensor from a noisy data tensor based on the least squares estimation with constant block approximation.
* `simulation`: Generate a symmetric tensor observation from the smooth signal tensor, Gaussian noise tensor, and permutation. Users can select one of 5 different smooth signal tensors generated from functions specified in Table 2.
* `simulation_bin`: Generate a symmetric binary tensor from the probability tensor and permutation. Users can select
one of 5 different smooth probability tensor generated from functions specified in Table 2.
* `simulation_asym`: Generate a non-symmetric tensor observation from the smooth signal tensor, Gaussian noise tensor, and permutation. Users can select one of 5 different smooth signal tensors generated from functions specified in Table S1.

`functions_asym.R`: assymetric version of `function_blk.R`

`tensor_visualization.R`: 3d tensor visualization function for Figure 5 and Figure S2

`Figure3.R`: R script for Figure 3 and Figure S3 in the main paper. Whole outputs of the simulation are in Data_Figure3.zip. It uses `functions_blk.R`

`Figure4.R`: R script for Figure 4 and Figure S1 in the main paper. Whole outputs of the simulation are in Data_Figure4.zip. It uses `functions_blk.R`


`Figure5.R`: R script for Figure 5 and Figure S2 in the main paper. It uses 'functions_blk.R' and `tensor_visualization.R`

`Table2.R`: R script for Table 2 computing the numerical CP and Tucker rank. It uses `functions_blk.R`

`Table3.R`: R script for Table 3 and Table S2. The simulation outputs are in Data_Table3.zip. It uses `functions_asym.R`

`revision.R`: R script for Figure R2 presented in revision`

Data_Figure3.zip: Simulation ouputs from the script `Figure3.R`

Data_Figure4.zip: Simulation ouputs from the script `Figure4.R`

Data_Table3.zip: Simulation ouputs from the script `Table3.R`
