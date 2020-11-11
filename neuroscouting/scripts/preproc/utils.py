from __future__ import division
import numpy as np
import nibabel as nib
import sys
import os


def load_nifti3D(filename):
    """ This function loads a 3D nifti file, converts it to a 3D array, and outputs the array matrix. Note that if you attempt to load a 4D nifti file this function will not work. """

    out=nib.load(filename)
    out_array=np.array(out.get_data())
    if len(out_array.shape)==4:
        # Since should be a 3D file, check that 4th dimension is only a single volume and if so, remove
        if out_array.shape[3]==1:
            out_array3D=out_array[:,:,:,0]
        else:
            sys.exit('4D file inputted, 3D file expected--CHECK!')
    elif len(out_array.shape)==3:
        # Already a 3D file
        out_array3D=out_array
    else:
        sys.exit('%dD file inputted, 3D file expected--CHECK!' %(len(out_array.shape)))

    return out_array3D


def ROIcheck(atlas_data,brainmask_data,propvoxels):
    """ This function masks an atlas file with a brainmask and loops through all atlas ROIs to determine whether proportion of voxels within the
 brainmask are greater than propvoxels. Prints to screen and outputs a list of invalid ROIs that are not sufficiently within the brainmask."""
    
    # Check to make sure mask and atlas dimensions match, otherwise exit out of this function
    if atlas_data.shape[0]==brainmask_data.shape[0] and atlas_data.shape[1]==brainmask_data.shape[1] and atlas_data.shape[2]==brainmask_data.shape[2]:

        # Determine number of ROIS. Masks are coded so each ROI has a numerical value, beginning with 1--NOTE THIS IS NOT NEEDED IF YOU HAVE ONLY ONE ROI IN A FILE
        nrois= int(len(np.unique(atlas_data))-1) #int(np.unique(atlas_data))

        # Mask atlas file with brainmask
        atlas_masked = np.multiply(atlas_data,brainmask_data)#--NOTE IF ALL ROIS ARE SEPARATE FILES, WILL NEED TO LOOP THROUGH THIS COMMAND SEPARATELY FOR EACH ONE HERE AND INCORPORATING THE BELOW, AS OPPOSED TO LOOPING THROUGH ROIS BELOW

        badrois=np.zeros(nrois)
        badcount=0
        badroilist=[]

        # Compare proportion of full atlas voxels with masked atlas voxels for each ROI
        if nrois == 116: # Only execute if working on the AAL atlas
            roiArr = np.unique(atlas_data[np.where(atlas_data>9000)])#[1:]
            for roiIDX, roi in enumerate(roiArr):
                if np.true_divide(np.where(atlas_masked==int(roi))[0].shape[0],np.where(atlas_data==int(roi))[0].shape[0])<propvoxels:# --PROPVOXELS I HAVE USUALLY DEFINED AS 0.25
                    badrois[roiIDX]=1 
                    badcount=badcount+1
                    badroilist.append(roiIDX) # Since this is the ROI value of 9001-9170
                else:
                    badrois[roiIDX]=0 
        else:
            for roi in range(1,nrois+1): # Add 1 so begins at value of 1 and not 0--THIS AND THE BADROIS DEFINITIONS WOULD BE OUTSIDE OF THE ATLAS_MASKED COMMAND IF EACH ROI IS A SEPARATE FILE
                if np.true_divide(np.where(atlas_masked==roi)[0].shape[0],np.where(atlas_data==roi)[0].shape[0])<propvoxels:# --PROPVOXELS I HAVE USUALLY DEFINED AS 0.25
                    badrois[roi-1]=1 # Subtract 1 because badrois begins indexing at 0, not 1
                    badcount=badcount+1
                    badroilist.append(roi) # Since this is the ROI value of 1-400, we will have to eventually subtract 1 from everything in this list if using python to make sure we index correctly
                else:
                    badrois[roi-1]=0 # Subtract 1 because badrois begins indexing at 0, not 1 

        print('\tMissing %d ROIs (fewer than %.2f voxels)' %(badcount,propvoxels))

    else:
        badrois=np.zeros(1)
        badcount=0
        badroilist=[]

    return badrois,badcount,badroilist


def get_sublist(subDir):
    """ Lists all directories in a main directories of all subjects to create a subject list. Note that this is BIDS format specific (assumes all subject directories begin "sub-xxx". Could be modified.
 Return: subject list """
    sublist = []
    for subnum in np.sort(os.listdir(subDir)): 
        if subnum[0:3] == 'sub': 
            if not 'html' in subnum: 
                sublist.append(subnum.split('-')[1]) 
    return sublist 

def get_seslist(subDir, sub):
    """ Lists all directories in a main directories of all subjects to create a subject list. Note that this is BIDS format specific (assumes all subject directories begin "sub-xxx". Could be modified.
 Return: subject list """
    seslist = []
    ses_topDir = os.path.join(subDir, 'sub-%s' %(sub))
    for sesnum in np.sort(os.listdir(ses_topDir)): 
        if sesnum[0:3] == 'ses': 
            seslist.append(sesnum.split('-')[1]) 
    return seslist 
