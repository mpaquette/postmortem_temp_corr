#!/usr/bin/env python3

import argparse
from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
import dipy.reconst.dti as dti
import nibabel as nib
import numpy as np
import pylab as pl


import warnings
warnings.simplefilter("ignore", RuntimeWarning)



description = """
Compute the volume-wise diffusivity correction factor from the DTI formulation [1].
Apply the correction to the log of normalized-signal and return it.
[1] Paquette et al., ISMRM 2024.
"""


def _build_args_parser():
    p = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('data', metavar='data', help='Path of the dwi nifti file.')
    p.add_argument('bval', metavar='bval', help='Path of the bval file.')
    p.add_argument('bvec', metavar='bvec', help='Path of the bvec file.')
    p.add_argument('index', metavar='index', help='Path to files with 1s in the position of volume to include')
    p.add_argument('out', metavar='output_image', help='Path of the output corrected data. (.nii.gz)')
    p.add_argument('out_ks', metavar='output_ks', help='Path of the diffusivity multiplier. (.txt)')
    p.add_argument('out_mult', metavar='output_mult', help='Path of the signal multiplier map. (.nii.gz)')
    p.add_argument('--mask', metavar='mask', help='Path of the brain mask for normalization.')
    
    return p



def main():
    parser = _build_args_parser()
    args = parser.parse_args()

    print('Loading data.')
    # load data
    img = nib.load(args.data)
    data = img.get_fdata()

    # load bval bvec
    bval, bvec = read_bvals_bvecs(args.bval, args.bvec)
    gtab = gradient_table(bval, bvec)
    # detect b0
    b0_th = 100.
    b0_index = np.where(bval < b0_th)[0]
    isnonb0 = bval > b0_th
    if bval.shape[0] != data.shape[3]:
        print('Data length different from bval length ({} vs {})'.format(data.shape[3], bval.shape[0]))
        return None 

    # Load index for DTI fitting.
    index = np.genfromtxt(args.index).astype(bool).ravel()
    if index.shape[0] != bval.shape[0]:
        print('index length different from bval length ({} vs {})'.format(index.shape[0], bval.shape[0]))
        return None 


    # check that index contains at least 1 b0
    if ~np.any(np.logical_and(~isnonb0, index)):
        print('Index must contain at least 1 b0.')
        return None


    if args.mask is None:
        mask = np.ones(data.shape[:3], dtype=bool)
    else:
        mask = nib.load(args.mask).get_fdata().astype(bool)

    totalVoxel = np.prod(mask.shape)
    voxelInMask = mask.sum()
    print('{} voxels out of {} inside mask ({:.1f}%)'.format(voxelInMask, totalVoxel, 100*voxelInMask/totalVoxel))






    data_untouched = data.copy()
    # step -1: Make sure data is clean
    data[np.isnan(data)] = 0
    data[np.isinf(data)] = 0
    data = np.clip(data, 0, np.inf)




    # Step 0: Compute spatial-mean data in mask
    data_spatial_mean = data[mask].mean(axis=0)
    #
    pl.figure()
    pl.subplot(2,2,1)
    pl.plot(np.arange(bval.shape[0])[isnonb0], data_spatial_mean[isnonb0], color='black', linewidth=2)
    pl.scatter(np.where(index), np.ones(index.sum())*(0.95*data_spatial_mean.min()), label='Included Volumes', color='blue', alpha=0.5)
    pl.legend()
    pl.xlabel('volume ordering')
    pl.title('Mean intensity (spatial, no b0s)')
    # pl.show()




    # Step 1: Compute diffusivity multiplier
    print('Computing calibration.') 
    ks = calibrate_many(data[mask], bvec, bval, index, current_k=None)
    #
    # pl.figure()
    pl.subplot(2,2,2)
    pl.plot(np.arange(bval.shape[0]), ks, color='black', linewidth=2)
    pl.scatter(np.where(index), np.ones(index.sum())*(0.99*ks.min()), label='Included Volumes', color='blue')
    pl.axhline(1.0, color='red', alpha=0.5, linestyle='dashed')
    pl.title('Diffusivity multipliers')
    pl.xlabel('volume ordering')
    # pl.show()
    #
    # Save Diffusivity multiplier coefs
    print('Saving Diffusivity multipliers')
    np.savetxt(args.out_ks, ks);




    # Step 2: Compute Mean b0 (from index)
    mean_b0 = data[..., np.logical_and(~isnonb0, index)].mean(axis=3)




    # Step 3: Normalize data and clean
    data_norm = data / mean_b0[..., None]
    data_norm[np.isnan(data_norm)] = 0
    data_norm[np.isinf(data_norm)] = 0
    data_norm = np.clip(data_norm, 0, 1)




    # Step 4: Compute Corrected data
    print('Computing corrected data')
    data_norm_corr = np.exp(np.log(data_norm)/ks)
    data_corr = mean_b0[..., None] * data_norm_corr
    #
    print('Saving Corrected data')
    nib.nifti1.Nifti1Image(data_corr.astype(np.float32), img.affine).to_filename(args.out)




    # Step 5: Compare spatial-mean data in mask
    data_corr_spatial_mean = data_corr[mask].mean(axis=0)
    #
    # pl.figure()
    pl.subplot(2,2,3)
    pl.plot(np.arange(bval.shape[0])[isnonb0], data_spatial_mean[isnonb0], color='black', linewidth=2, label='before', alpha=0.75)
    pl.plot(np.arange(bval.shape[0])[isnonb0], data_corr_spatial_mean[isnonb0], color='red', linewidth=2, label='after', alpha=0.75)
    pl.legend()
    pl.xlabel('volume ordering')
    pl.title('Mean intensity (spatial, no b0s)')
    # pl.show()




    # Step 6: Compute Voxel-wise signal multiplier
    print('Computing voxel-wise signal multipliers')
    signal_mult = data_corr / data
    signal_mult[np.isnan(signal_mult)] = 1
    signal_mult[np.isinf(signal_mult)] = 1
    signal_mult = np.clip(signal_mult, 0, np.inf)
    print('Saving voxel-wise signal multipliers')
    nib.nifti1.Nifti1Image(signal_mult, img.affine).to_filename(args.out_mult)




    # Step 7: Compute mean mult in mask
    signal_mult_spatial_mean = signal_mult[mask].mean(axis=0)
    signal_mult_spatial_median = np.median(signal_mult[mask], axis=0)
    #
    # pl.figure()
    pl.subplot(2,2,4)
    pl.plot(np.arange(bval.shape[0]), signal_mult_spatial_mean, color='black', linewidth=2, label='mean')
    pl.plot(np.arange(bval.shape[0]), signal_mult_spatial_median, color='blue', linewidth=2, label='median')
    pl.legend()
    pl.xlabel('volume ordering')
    pl.title('Signal multiplier (spatial)')
    pl.show(block=False)
    pl.pause(30) # close after 30sec if no interaction
    pl.close()



def calibrate_many(data, bvecs, bvals, keep_idx, current_k=None):
    if current_k is None:
        current_k = np.ones(bvals.shape[0])
    # apply ks on the bvals to trick the fit
    bvals_mod = bvals * current_k
    # 
    gtab = gradient_table(bvals=bvals, bvecs=bvecs) 
    gtab_low_mod = gradient_table(bvals=bvals_mod[keep_idx], bvecs=bvecs[keep_idx])
    #
    dti_model_low_mod = dti.TensorModel(gtab_low_mod, fit_method='WLS', return_S0_hat=True)
    #
    fit_low_mod = dti_model_low_mod.fit(data[:, keep_idx]) # this estimates D_ss from heated signal because its using bval_mod
    #
    S0_hat = fit_low_mod.S0_hat
    S_pred = fit_low_mod.predict(gtab) # this is estimates S_ss because its using estimated D_ss with normal bval
    S_pred_norm = S_pred / S0_hat[:, None]
    #
    S_pred_norm[np.isnan(S_pred_norm)] = 0
    S_pred_norm[np.isinf(S_pred_norm)] = 0
    #
    pred_adc = -np.log(S_pred_norm) / bvals # this is estimating adc_ss because S_pred is estimating S_ss because we divide by normal bval
    #
    pred_adc[np.isinf(pred_adc)] = 0
    pred_adc[np.isnan(pred_adc)] = 0
    #
    ks = 1 - (np.log(data/S_pred)/(bvals*pred_adc)) # this estimate full ks because pred and adc are estimating ss quantities  
    #
    ks[np.isnan(ks)] = 1
    ks[np.isinf(ks)] = 1
    #
    return np.median(ks, axis=0)





if __name__ == "__main__":
    main()






