# Workflow
### Step 1. Generate data
- expression dataset
- a continuous modulator
- 
A varying coefficient model is considered, where the regression coefficients vary across samples according to the modulator.  
Based on this model, the response is generated with additive Gaussian noise.

### Step 2. Build a hyperparameter grid
A grid is created over:
- `lambda`
- `alpha`
- `h`
For each target sample, every hyperparameter combination is evaluated using `ssGCV()`.

### Step 3. Select sample-specific hyperparameters and fit a kernel-based elastic net model
For each target sample, the hyperparameter combination with the smallest `ssGCV` value is selected.  
Using the selected hyperparameters, a kernel-based elastic net model is then fitted.

### Step 4. Interpret outputs
Outputs are: 
- `opt`_mx: selected `lambda`, `alpha`, and `h` for each target sample.
- `BETA`: sample-specific coefficient estimates.

---

# `ssGCV()` function
### Arguments
- X: predictor matrix
- y: response vector
- mod_vec: modulator values
- target: index of the target sample
- lambda: elastic net regularization parameter
- alpha: elastic net mixing parameter
- bandwidth: Gaussian kernel bandwidth

### Returns
- sample-specific GCV score for the target sample under the given hyperparameter setting
