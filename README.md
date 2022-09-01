# Model-Resolution-based-Deconvolution-for-FD-Photoacoustic-Reconstruction
Model resolution based deconvolution (SALSA) improves frequency domain photo acoustic reconstruction. The full paper for reference can be found in https://doi.org/10.1121/10.0013829

## Step 1 (Reconstruction_penalty_function_step1.m)
### Parametres section
Input the parameters in _Paramters_ section. Range of frequencies selected - (100kHz - 2.5Mhz) for _FUNDUS and Derenzo Phantom_ and _(5kHz - 150kHz)_ for _Breast phantom_. Set the required SNR and FOV. For split/ discontinuous frequency ranges, set _fh_end, fh_start, wnh_ parameters.

### Calculate W Matrix
Calculate Weight matrix (forward modelling) using function _Calculate_W_matrix.m_ function. Regularization paramter for _Derenzo phantom and FUNDUS is 100_ and _2e11 for Breast phantom_.

### Calculate pressure matrix
Calculte _P = WX_ to get pressure matrix (Complete forward model).

### Inversion model using different penalty functions
Reconstruction is performed with a penalty function multiplied to regularization paramter.

# Step 2 (Model_resolution_step2.m)
Set the same paramters as in step 1.
All the deconvolutions were based on **SALSA** approach and were done for 100 iterations for each of the penalized reconstructions.
