import argparse
import os
import sys
import glob
import pandas
import numpy
from subprocess import call
import collections

parser = argparse.ArgumentParser(prog="ROIExtractor")
parser.add_argument("input_BIDS", help="target BIDS dataset that has derivatives and fmriprep preprocessed files present")
parser.add_argument("-task_name", help="task name to analyze.", default = 0)
parser.add_argument("-targetsuffix", help="target image file suffix. Defaults to Asym_preproc.nii.gz", default = "Asym_preproc.nii.gz")
parser.add_argument("-session", help="session to analyze. (0 analyzes all sessions)", default = 0)
parser.add_argument("-subject", help= "subject to analyze.")
parser.add_argument("-input_fmriprep", help= "location of the fmriprep output folder.")
group = parser.add_mutually_exclusive_group()
group.add_argument("-roimask", help= "nifti file that contains the ROI mask as integers")

parser.add_argument("-brainmask", help= "subject brain mask suffix. Defaults to brainmask.nii.gz. Highly recommended to ensure voxels outside of brain mask are not included in average", default = "brainmask.nii.gz")
parser.add_argument("-atlasName", help = "Name of the output ROI time series (default: rois)", default = "rois")
parser.add_argument("--nodatacheck", help="if option is included, will not inform user how many datasets are about to be processed, and will not require final confirmation", action = "store_true")
parser.add_argument("--noregcheck", help="if option is included, will not inform user what regression/method is used, and will not require final confirmation", action = "store_true")
parser.add_argument("--overwrite", help="allow for overwriting regressor matrices and output nii.gz", action="store_true")
parser.add_argument("-q","--quiet", help="do NOT print activity to console, activity will still be recorded in log file", action = "store_true")

args = parser.parse_args()

targetdirs = args.input_BIDS
fmriprepdirs = args.input_fmriprep

print('TARGET DIRECTORY: ', targetdirs)
if not os.path.isdir(targetdirs):
    print("fmriprep directory is not present in the BIDS/derivatives directory. Have you run fmriprep? Ending function...")
    sys.exit()

files = glob.glob(targetdirs+"/**/*.nii.gz", recursive = True)

targets = [i for i in files if "bold" in i]
print(args.targetsuffix)
targets = numpy.sort([i for i in targets if args.targetsuffix in i])
if args.task_name is not 0:
    targets = [i for i in targets if "task-"+args.task_name +"_" in i]
if args.session is not 0:
    targets = [i for i in targets if "ses-" + args.session + "_" in i]
if args.subject is not None:
    targets = [i for i in targets if "sub-" + args.subject + "_" in i]
print(targets)
subs = [i.split("/") for i in targets]
subList = []
for file in subs:
    subtarget = [i for i in file if "sub-" in i]
    subList.extend([subtarget[0]])
uniqueSubs = set(subList)

counter = collections.Counter(subList).values()
minNum = min(counter)
maxNum = max(counter)
dataDesc = "This data set contains " + str(len(uniqueSubs)) + " subjects with at max " +          str(maxNum) +" and at minimum " + str(minNum) +" functional scans. Does this look correct?"

if not args.nodatacheck:
    print(dataDesc)
    while True:
        userinput = input("y/n: ")
        if userinput == "n":
            print("Ending this script. Please check your dataset and try again.")
            sys.exit()
        if userinput != "y" and userinput !="n":
            continue
        else:
            break

for target in targets:
        output = []
        bits = target.split("_")
        for bidx in bits:
            if 'sub-' in bidx:
                b2=bidx.split("/")
                sub=b2[len(b2)-1]
            if 'ses-' in bidx and len(bidx) == 7:
                ses = bidx#.split('-')[1]
            if 'run-' in bidx and len(bidx) == 7:
                run = bidx#.split('-')[1]
            if '.nii.gz' in bidx:
                BPF = bidx[:-len('.nii.gz')]
        print(sub)
        print(fmriprepdirs + '/' + sub + '/func' + '/*' + run + '*' + args.brainmask)
        sub_brainmask = glob.glob(fmriprepdirs + '/' + sub + '/func' + '/*' + run + '*' + args.brainmask)
        print(sub_brainmask)

        comTarget = bits.index("bold.nii.gz")+4
        bits = bits[0:(comTarget)]
        targetRegs = '_'.join(bits)
        targetOutput = targetRegs+'_'+BPF+"_"+args.atlasName+".tsv"
        print('ROI MASK: ', args.roimask)
        print('BRAIN MASK: ', sub_brainmask[0])
        fslCommand = "fslmeants -i " + target + " -m " + sub_brainmask[0] + " --label=" + args.roimask + " -o " +targetOutput
        print(fslCommand, '\n')
        call(fslCommand, shell =True)
