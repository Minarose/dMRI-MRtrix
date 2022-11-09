#!/bin/bash

#display_usage() {
  #echo "$(basename $0) [b1000 ni/nii.gz] [b2000 ni/nii.gz] [b3000 ni/nii.gz] "
  #echo "This script uses MRtrix to analyze diffusion data. It requires 3 arguments:
    #1) b1000 ni/nii.gz;
    #2) b2000 ni/nii.gz;
    #3) b3000 ni/nii.gz
    #4)"
  #}

  #if [ $# -le 2 ]
  #then
    #display_usage
    #exit 1
  #fi

#b100=$1
#b200=$2
#b300=$3

########################### STEP 1 ###################################
#             Convert data to .mif format and denoise                        #
######################################################################

#Warning PAR/REC files must be converted to NIFTI

#cd to directory with NIFTI files
cd /home/minarose/Desktop/BIDS_ATTEMPT/NIFTI/DWI

#Convert b1000/2000/3000/0 NIFTI to mif files using mrtrix
mrconvert sub-01_b1000.nii.gz sub-01_b1000.mif -fslgrad sub-01_b1000.bvec sub-01_b1000.bval #first variable
mrconvert sub-01_b2000.nii.gz sub-01_b2000.mif -fslgrad sub-01_b2000.bvec sub-01_b2000.bval #second variable
mrconvert sub-01_b3000.nii.gz sub-01_b3000.mif -fslgrad sub-01_b3000.bvec sub-01_b3000.bval #third variable
mrconvert sub-01_b0.nii.gz sub-01_b0.mif  #fourth variable

#concatenate all b mif images to one dwi mif file
mrcat sub-01_b1000.mif sub-01_b2000.mif sub-01_b3000.mif raw_dwi.mif

#denoise
dwidenoise raw_dwi.mif dwi_den.mif -noise noise.mif

#calculate difference between raw and denoised image
mrcalc raw_dwi.mif dwi_den.mif –subtract residual.mif 
##mrview noise.mif residual.mif #outputs a .mif noise image

#Gibbs denoising (using mrdegibbs) to check your diffusion data for ringing
mrdegibbs dwi_den.mif dwi_den_unr.mif –axes 0,1

#calculate difference between denoised image and unringed image
mrcalc dwi_den.mif dwi_den_unr.mif –subtract residualUnringed.mif 
##mrview dwi_den_unr.mif residualUnringed.mif

#################################STEP 2################################
#				Preprocessing					#
#######################################################################

###################EPI Distortion Correction############################
#for PA direction, extract all b0s from dwi_den_unr.mif and calculate mean
dwiextract dwi_den_unr.mif - -bzero | mrmath – mean mean_b0_AP.mif –axis 3 

#calculate mean for PA direction
mrconvert sub-01_b0.mif - | mrmath – mean mean_b0_PA.mif –axis 3

#concatenate two mean b0 images together: mrcat mean_b0_AP.mif mean_b0_PA.mif –axis 3 b0_pair.mif
mrcat mean_b0_AP.mif mean_b0_PA.mif –axis 3 b0_pair.mif

###################TOPUP/EDDY#############################################
dwifslpreproc dwi_den_unr.mif dwi_den_unr_preproc.mif -nocleanup -pe_dir AP -rpe_pair -se_epi b0_pair.mif -eddy_options " --slm=linear --data_is_shelled"
#might take 2-4 hours (13:21-14:34)
##################BIAS FIELD CORRECTION###################################

dwibiascorrect ants dwi_den_unr_preproc.mif dwi_den_unr_preproc_unbiased.mif -bias bias.mif
#if ants unavailbl use fsl
#dwibiascorrect fsl dwi_den_unr_preproc.mif dwi_den_unr_preproc_unbiased.mif –bias bias.mif 

#check results
mrview dwi_den_unr_preproc.mif & mrview dwi_den_unr_preproc_unbiased.mif & mrview bias.mif –colourmap 2

##################BRAIN MASK ESTIMATION #######################################
dwi2mask dwi_den_unr_preproc_unbiased.mif mask_den_unr_preproc_unb.mif 

#a. if mask is overly-inclusive, refine mask using maskfilter
#maskfilter -scale 5 mask_den_unr_preproc_unb.mif clean mask_den_unr_preproc_unb_filt.mif -nthreads 4 

#note: scale needs to be adjusted for each subject 
# check results: 
#mrview dwi_den_unr_preproc_unbiased.mif -overlay.load mask_den_unr_preproc_unb_filt.mif

########################### STEP 3 ###################################
#             Basis function for each tissue type                    #
######################################################################
# Create a basis function from the subject's DWI data. The "dhollander" function is best used for multi-shell acquisitions; it will estimate different basis functions for each tissue type. For single-shell acquisition, use the "tournier" function instead
dwi2response dhollander dwi_den_unr_preproc_unbiased.mif wm.txt gm.txt csf.txt -voxels voxels.mif -nthreads 4 
#view response functions for WM, GM, CSF on shview wm.txt; shview gm.txt shview csf.txt?

#Performs multishell-multitissue constrained spherical deconvolution, using the basis functions estimated above
dwi2fod msmt_csd dwi_den_unr_preproc_unbiased.mif –mask mask_den_unr_preproc_unb.mif wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif -nthreads 4
# Creates an image of the fiber orientation densities overlaid onto the estimated tissues (Blue=WM; Green=GM; Red=CSF)
# You should see FOD's mostly within the white matter. These can be viewed later with the command "mrview vf.mif -odf.load_sh wmfod.mif"
mrconvert –coord 3 0 wmfod.mif - | mrcat csffod.mif gmfod.mif – vf.mif 
mrview vf.mif -odf.load_sh wmfod.mif

# Now normalize the FODs to enable comparison between subjects
mtnormalise wmfod.mif wmfod_norm.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif –mask mask_den_unr_preproc_unb.mif


########################### STEP 3 ###################################
#            Create a GM/WM boundary for seed analysis               #
######################################################################

#T1 to DW coregistration
#preprocess T1 image
bet sub-01_T1.nii.gz T1_bet.nii.gz
5ttgen fsl -premasked T1_bet.nii.gz 5tt_nocoreg.mif -nthreads 4 
#*** error: qform and sform are inconsistent in NTfTI image T1.nii - using sform 
#b. coregistration: 
dwiextract dwi_den_unr_preproc_unbiased.mif - -bzero | mrmath – mean mean_b0_preprocessed.mif –axis 3 mrconvert mean_b0_preprocessed.mif mean_b0_preprocessed.nii 

flirt -in mean_b0_preprocessed.nii -ref T1_bet.nii.gz -dof 6 -omat diff2struct_fsl.mat 
transformconvert diff2struct_fsl.mat mean_b0_preprocessed.nii T1_bet.nii.gz flirt_import diff2struct_mrtrix.txt 
mrconvert T1_bet.nii.gz T1_raw.mif 
mrtransform T1_raw.mif –linear diff2struct_mrtrix.txt –inverse T1_coreg.mif 
mrtransform 5tt_nocoreg.mif –linear diff2struct_mrtrix.txt –inverse 5tt_coreg.mif 
#c. check results: 
mrview dwi_den_unr_preproc_unbiased.mif –overlay.load T1_raw.mif –overlay.colourmap 2 -overlay.load T1_coreg.mif –overlay.colourmap 1 
#d. check 5tt_coreg results: 
mrview dwi_den_unr_preproc_unbiased.mif -overlay.load 5tt_coreg.mif

########################### STEP 4 ###################################
#                 Run the streamline analysis                        #
######################################################################

#Prepare streamline mask
5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif

mrview dwi_den_unr_preproc_unbiased.mif -overlay.load gmwmSeed_coreg.mif

# Create streamlines
# Note that the "right" number of streamlines is still up for debate. Last I read from the MRtrix documentation,
# They recommend about 100 million tracks. Here I use 10 million, if only to save time. Read their papers and then make a decision
tckgen –act 5tt_coreg.mif –backtrack –seed_gmwmi gmwmSeed_coreg.mif –select 1000000 wmfod_norm.mif tracts_1million.tck -nthreads 4

#note: when viewing streamlines reduce number to 100k: 
tckedit tracts_1million.tck –number 100k tracts_100k.tck 
#a. check results: 
mrview dwi_den_unr_preproc_unbiased.mif –tractography.load tracts_100k.tck 
#note: change slab thickness to 5/10 when viewing tracts


# Reduce the number of streamlines with tcksift
tcksift -act 5tt_coreg.mif -term_number 500000 tracts_1million.tck wmfod_norm.mif sift_500k.tck 
#a. edit to 100k tracts: 
tckedit sift_500k.tck -number 100k tracts_sifted_100k.tck 
#b. check results: 
mrview dwi_den_unr_preproc_unbiased.mif -tractography.load tracts_sifted_100k.tck

########################################################

# These commands are for quality-checking your diffusion data
#most already done earlier sort them later

### Quality checks for Step 2 ###

# Views the voxels used for FOD estimation
echo "Now viewing the voxels used for FOD estimation (Blue=WM; Green=GM; Red=CSF)"
mrview dwi_den_unr_preproc_unbiased.mif -overlay.load voxels.mif

# Views the response functions for each tissue type. The WM function should flatten out at higher b-values, while the other tissues should remain spherical
echo "Now viewing response function for white matter (press right arrow key to view response function for different shells)"
shview wm.txt
echo "Now viewing response function for grey matter"
shview gm.txt
echo "Now viewing response function for CSF"
shview csf.txt

# Views the FODs overlaid on the tissue types (Blue=WM; Green=GM; Red=CSF)
mrview vf.mif -odf.load_sh wmfod.mif

### Quality checks for Step 3 ###

# Check alignment of the 5 tissue types before and after alignment (new alignment in red, old alignment in blue)
mrview dwi_den_unr_preproc_unbiased.mif -overlay.load 5tt_nocoreg.mif -overlay.colourmap 2 -overlay.load 5tt_coreg.mif -overlay.colourmap 1

# Check the seed region (should match up along the GM/WM boundary)
mrview dwi_den_unr_preproc_unbiased.mif -overlay.load gmwmSeed_coreg.mif


### Quality checks for Step 4 ###

# View the tracks in mrview
mrview dwi_den_unr_preproc_unbiased.mif -tractography.load tracts_sifted_100k.tck 

### Quality checks for Step 3 ###

# Check alignment of the 5 tissue types before and after alignment (new alignment in red, old alignment in blue)
mrview dwi_den_unr_preproc_unbiased.mif -overlay.load 5tt_nocoreg.mif -overlay.colourmap 2 -overlay.load 5tt_coreg.mif -overlay.colourmap 1

# Check the seed region (should match up along the GM/WM boundary)
mrview dwi_den_unr_preproc_unbiased.mif -overlay.load gmwmSeed_coreg.mif


### Quality checks for Step 4 ###

# View the tracks in mrview
mrview dwi_den_unr_preproc_unbiased.mif -tractography.load tracts_sifted_100k.tck


###################################################################
#                       Preprocess T1		                  #
###################################################################

mrconvert T1_raw.mif T1_raw.nii.gz 

recon-all -s sub-1 -i T1_raw.nii.gz -all

########################### STEP 5 ###################################
#                 Prepare atlas for structural connectivity          #
######################################################################
#Download annotation file from here: https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446
#copy files into $SUBJECTS_DIR/fsaverage/label/ 

#Map the annotation files of the HCP MMP 1.0 atlas from fsaverage to you subject. Remember to do that for both hemispheres:

mri_surf2surf --srcsubject fsaverage --trgsubject sub-1 --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.HCP-MMP1.annot --tval $SUBJECTS_DIR/sub-1/label/lh.HCP-MMP1.annot

mri_surf2surf --srcsubject fsaverage --trgsubject sub-1 --hemi rh --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.HCP-MMP1.annot --tval $SUBJECTS_DIR/sub-1/label/rh.HCP-MMP1.annot 

#download files of two colour-lookup tables from supplementary files in BATMAN tutorial file (hcpmmp1_original.txt & hcpmmp1_ordered.txt) to fsaverage/label directory too? from: https://osf.io/fkyht/
 

#Map the HCP MMP 1.0 annotations onto the volumetric image and add (FreeSurfer-specific) subcortical segmentation.: 
mri_aparc2aseg --old-ribbon --s sub-1 --annot HCP-MMP1 --o HCP-MMP1.mgz

#Convert the resulting file to .mif format (use datatype uint32, which is liked best by MRtrix): 
mrconvert -datatype uint32 HCP-MMP1.mgz HCP-MMP1.mif 

#Replace the random integers of the hcpmmp1.mif file with integers that start at 1 and increase by 1:
labelconvert HCP-MMP1.mif $SUBJECTS_DIR/fsaverage/label/hcpmmp1_original.txt $SUBJECTS_DIR/fsaverage/label/hcpmmp1_ordered.txt hcpmmp1_parcels_nocoreg.mif 


#Register the ordered atlas-based volumetric parcellation to diffusion space: 
mrtransform hcpmmp1_parcels_nocoreg.mif –linear diff2struct_mrtrix.txt –inverse –datatype uint32 hcpmmp1_parcels_coreg.mif


########################### STEP 6 ###################################
#                 Matrix generation			             #
######################################################################

#create .csv file containing matrix: 
tck2connectome –symmetric –zero_diagonal -scale_invnodevol sift_500k.tck hcpmmp1_parcels_coreg.mif hcpmmp1.csv –out_assignment assignments_hcpmmp1.csv 

#b. visualize matrix using MATLAB: import data --> hcpmmp1.csv --> output type --> numeric matrix --> import selection --> import data --> imagesc (hcpmmp1) --> reducedlimit=[0 0.1] --> imagesc (hcpmmp1, reducedlimit) --> axis square

#importdata(filename), reducedlimit = [0 0.1], imagesc(csv, reducedlimit) can probably be done from terminal with matlab comands right? or do we need to open matlab



