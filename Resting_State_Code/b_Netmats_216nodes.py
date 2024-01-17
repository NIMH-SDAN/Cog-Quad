#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by Anderson Winkler 02/25/2019
Modified by Julia Linke 05/01/2020

"""
#the script assumes that there is a directory called lists in the rootdir
#that contains a list with the participant IDs 

import nibabel
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import imageio

rootdir = os.path.join('/','data','SMDN','RestCogQuad')
atlasdir = os.path.join('/','data','SMDN','CHANGE2')
partial = False
remove_vols = 0
fd_thr = .5

def loadGifti(fname, NonSteadyState=0, icres=7): # ============================
    gii = nibabel.load(fname)
    gii_data = [d.data[:,None] for d in gii.darrays]
    gii_data = np.concatenate(gii_data, axis=1).T
    nV = 4**icres*10 + 2
    gii_data = gii_data[:,0:nV]
    return gii_data[NonSteadyState:,:]

def loadNifti(fname, NonSteadyState=0): # =====================================
    n = nibabel.load(fname)
    naff   = n.affine
    img4d  = n.get_fdata();
    imgsiz = img4d.shape
    if len(imgsiz) == 4:
        img4d  = img4d[:,:,:,NonSteadyState:]
        imgsiz = img4d.shape
        img2d  = np.reshape(img4d, (np.prod(imgsiz[0:3]), imgsiz[-1]), order='F').T
    else:
        img2d  = np.reshape(img4d, (np.prod(imgsiz[0:3]), 1), order='F').T
    return img2d, imgsiz, naff

def unWrap(netmat, side='lower'): # ===========================================
    if side.lower() == 'lower':
        uw  = netmat[np.tril_indices(netmat.shape[0], k=-1, m=netmat.shape[1])]
    elif side.lower() == 'upper':
        uw  = netmat[np.triu_indices(netmat.shape[0], k= 1, m=netmat.shape[1])]
    elif side.lower() == 'both':
        uwl = netmat[np.tril_indices(netmat.shape[0], k=-1, m=netmat.shape[1])]
        uwu = netmat[np.triu_indices(netmat.shape[0], k= 1, m=netmat.shape[1])]
        uw  = np.concatenate((uwl, uwu))
    return uw

# ======== [ Main ] ===========================================================
if __name__ == '__main__':
    
    #We could run the next steps for all participants in the fmriprep directory like this:
    #listSubjects = [d for d in os.listdir(os.path.join(rootdir,'derivatives','fmriprep')) if (os.path.isdir(os.path.join(rootdir,'derivatives','fmriprep',d)) and d[0:4] == 'sub-')]
    #However, we run it for only participants in our list
    listSubjects = open(os.path.join(rootdir,'lists/FinalSample2.txt'),'r')
    listSubjects = listSubjects.read().split()
    listSubjects = [x.replace('s','sub-s') for x in listSubjects]
    listSubjects.sort()
    fmriprepdir  = os.path.join(rootdir, 'derivatives/fmriprep/')
    confounddir  = os.path.join(rootdir,'derivatives/confounds')
    netmatsdir   = os.path.join(rootdir,'derivatives/netmats')
    
    # List of subcortical structures of interest
    # 26/58 = Accumbens
    # 18/54 = Amygdala
    # 11/50 = Caudate
    # 17/53 = Hippocampus
    # 13/52 = Pallidum
    # 12/51 = Putamen
    # 10/49 = Thalamus
    # 28/60 = Ventral DC
    aseg_list = {
            'L': [26, 18, 11, 17, 13, 12, 10, 28],
            'R': [58, 54, 50, 53, 52, 51, 49, 60]}
    
    # Load the parcellation "annot" files in fsaverage space
    annot = {}
    ctab  = {}
    names = {}
    [annot['L'], ctab['L'], names['L']] = nibabel.freesurfer.io.read_annot(os.path.join(atlasdir ,'atlas', 'Schaefer2018', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label','lh.Schaefer2018_200Parcels_17Networks_order.annot'))
    [annot['R'], ctab['R'], names['R']] = nibabel.freesurfer.io.read_annot(os.path.join(atlasdir ,'atlas', 'Schaefer2018', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label','rh.Schaefer2018_200Parcels_17Networks_order.annot'))
    
    # For each subject
    # Maybe we need to select a subset of participants
    # for sidx, subj in enumerate(listSubjects[40:]):
    # Right now we just take all participants
    for sidx, subj in enumerate(listSubjects):
        print(sidx, subj)
        
        # Load the confounds; unfortunately, there are 2 different delimiters
        # in the confounds file > need to accommodate that
        Rest_confounds = pd.read_csv(os.path.join(confounddir, subj,'{}_ses-1_task-rest_AllConfounds.tsv'.format(subj)), header=None, sep='\t| |[|]',engine='python')
        
        Rest_confounds = Rest_confounds.values
        
        # For each hemisphere
        Rest_surf_func = {}
        Rest_surf_parc = {}
        
        for hemi in ['L', 'R']:
                    
            # Load functional data in fsaverage space (surface)
            Rest_surf_func[hemi] = loadGifti(os.path.join(fmriprepdir, subj, 'ses-1', 'func', '{}_ses-1_task-rest_space-fsaverage_hemi-{}_bold.func.gii'.format(subj, hemi)))
            
            # Regress out confounds
            Rest_b = np.linalg.lstsq(Rest_confounds, Rest_surf_func[hemi], rcond=None)[0]
            Rest_surf_func[hemi] = Rest_surf_func[hemi] - np.matmul(Rest_confounds, Rest_b)
            
            # For each cortical parcel, extract the average timecourse
            U = np.unique(annot[hemi])
            Rest_surf_parc[hemi] = np.zeros((Rest_surf_func[hemi].shape[0], U.shape[0]))
            
            for parc in U:
                Rest_surf_parc[hemi][:,parc] = np.mean(Rest_surf_func[hemi][:,annot[hemi] == parc], axis=1)
                     
        # Load the subcortical segmentation (FS "aseg" files) in MNI space
        [Rest_aseg, Rest_aseg_siz, Rest_aseg_aff]   = loadNifti(os.path.join(fmriprepdir, subj,'ses-1','func', '{}_ses-1_task-rest_space-T1w_desc-aseg_dseg.nii.gz'.format(subj)))
        
        # Load their functional data in MNI space (volume)
        [Rest_vol_func, Rest_vol_siz, Rest_vol_aff] = loadNifti(os.path.join(fmriprepdir, subj,'ses-1','func', '{}_ses-1_task-rest_space-T1w_desc-preproc_bold.nii.gz'.format(subj))) 
        
        # Regress out confounds
        Rest_b = np.linalg.lstsq(Rest_confounds, Rest_vol_func, rcond=None)[0]
        Rest_vol_func = Rest_vol_func - np.matmul(Rest_confounds, Rest_b)      
        
        # For each subcortical parcel, extract the average timecourse
        Rest_vol_parc = {}
        for hemi in ['L', 'R']:
            Rest_vol_parc[hemi] = np.zeros((Rest_vol_func.shape[0], len(aseg_list[hemi])))
            for pidx, parc in enumerate(aseg_list[hemi]):
                Rest_vol_parc[hemi][:,pidx] = np.mean(Rest_vol_func[:,np.squeeze(Rest_aseg == parc)], axis=1)
     

        # Merge cortical and subcortical timecourses of both hemispheres
        Rest_all_parc = np.concatenate((Rest_surf_parc['L'][:,1:], Rest_surf_parc['R'][:,1:], Rest_vol_parc['L'], Rest_vol_parc['R']), axis=1)
        
        # Compute the "netmat" between the regions (cortical and subcortical),
        # and unwrap it      
        Rest_rmat = np.corrcoef(Rest_all_parc, rowvar=False);
        if partial == True:
            Rest_rinv = np.linalg.pinv(Rest_rmat)
            Rest_diag = np.diagonal(Rest_rinv)[:,None]
            Rest_rmat = np.multiply(-Rest_rinv, np.power(np.multiply(Rest_diag, Rest_diag.T), -.5))
            
        # Save the wrapped matrices as csv file 
        if not os.path.exists(os.path.join(netmatsdir, subj)):
            os.makedirs(os.path.join(netmatsdir, subj))
        np.savetxt(os.path.join(netmatsdir, subj,'{}_ses-1_task-rest_netmat-{}_atlas-Schaefer2018-200P+17N_space-T1w.csv'.format(subj, 'partial' if partial else 'full')), Rest_rmat, delimiter=",")
        imageio.imwrite(os.path.join(netmatsdir, subj,'{}_ses-1_task-rest_netmat-{}_atlas-Schaefer2018-200P+17N_space-T1w.png'.format(subj, 'partial' if partial else 'full')), (255*(Rest_rmat+1)/2).astype('uint8'))
        plt.imshow(Rest_rmat)
        plt.show()
        
        
        
        