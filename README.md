# Correcting temperature related diffusivity drift for postmortem difusion MRI
## Michael Paquette, Cornelius Eichner, Christian Bock, and Alfred Anwander
### ISMRM2024 
#### Program 2423

This repo contains a short example of the presented method using invivo data with artificial diffusivity corruption.  
The dataset was created from a single b=1000 shell from subject PT001_ses-1_acq-1 from the [Pentera 3T public dataset](https://zenodo.org/records/2602049)  


The Ground truth data (1 b0 + 32 b1000) was corrupted to simulate temperature induced diffusitity drift.  

![Spherical means of the ground truth normalized data, the arificially corrupted data and the corrected data.](images/spherical_means.png)

We use the b0 and last 16 volumes ([see index file](corrupted_data/index.txt)) to estimate the steady-state diffusivities.  

```bash
dti_temp_correction.py corrupted_data/corrupted.nii.gz \
                       data/bval.txt \
                       data/bvec.txt \
                       corrupted_data/index.txt \
                       corrected_data/corrected.nii.gz \
                       corrected_data/estimated_coef.txt \
                       corrected_data/signal_multipliers.nii.gz \
                       --mask data/mask.nii.gz
```

![Ground truth vs estimated alpha coefficients.](images/alpha_coef_estimation.png)

![Mean squared error for corrupted and corrected data with respect to the ground truth.](images/signal_mse.png)
