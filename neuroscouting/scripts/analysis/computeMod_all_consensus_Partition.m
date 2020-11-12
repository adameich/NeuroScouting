% The steps are:
% - Run the Louvain community detection algorithm multiple times, saving the output on each iteration. Use [community_louvain.m]
% - Take the saved community affiliations from each run and turn them into an affiliation vector. Use [agreement.m]
% - Create a consensus partition, running it nReps times, with a specified threshold of agreement tau. Use [consensus_und.m]
% - Given the consensus partition, compute the modularity value a final time. Use [community_louvain.m] again I think


% 1. Compute modularity N times and save all output in an NNodes x N mod reps matrix
% 2. Create agreement matrix (a proportion, (0-N)/N, of times a node is in the same module)
% 3. Get consensus partition, running Nreps times, with a threshold of tau
% 4. Calculate modularity given the partition
% 5. Save output partitions, modularity, and number modules


clear all

testBug = 0;

topDir = '/Users/Eichenbaum/HWNI';

%% Set up path structures for MATLAB
addpath /Users/Eichenbaum/HWNI/software/BCT/2017_01_15_BCT;
addpath /Users/Eichenbaum/HWNI/software/npy-matlab-master;
addpath /Users/Eichenbaum/HWNI/Experiments/nScout/scripts/analysis/MatlabUtils;


%% Depath paths and variables
basedir = '/Users/Eichenbaum/HWNI/Experiments/nScout/';
scriptdir = [basedir 'scripts/'];
subdir = [basedir 'data/derivatives/denoise_out/'];
outdir = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/derivatives/group-level/';

atlas = 'brainnetome'; %'Schaefer400node' 'brainnetome'

if strcmp(atlas, 'Schaefer400node')
    missingROIsFName = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/analyses/missingROIs_list_thresh-0.25_Schaefer400node.npy';
    commonROIsFName = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/analyses/commonROIs_list_thresh-0.25_Schaefer400node.npy';
elseif strcmp(atlas, 'brainnetome')
    missingROIsFName = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/analyses/missingROIs_list_thresh-0.25_Brainnetome.npy';
    commonROIsFName = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/analyses/commonROIs_list_thresh-0.25_Brainnetome.npy';
end

missingROIs = readNPY(missingROIsFName);
commonROIs = readNPY(commonROIsFName);

%Parameters for consensus clustering (numbers taken from JRCs script in which she says "taken from my interpretation of Lancichinette & Fortunato, 2012")
Nmods = 1000; % Number of times optimal participations are generated
Nreps = 1000; % Number of times consensus cluster is partitioned
tau = 0.5; % Proportion of times a node must be in the same partition to be counted

% % Proportional threshold value (usually 0.05, 0.075, 0.10, 0.125, or 0.15)
% prop_thresh = 0.125;

subList = get_sublist(subdir);
% subList = subList(1:69); %% Only for Lab Meeting

if testBug;
    subList = {'sub-201'};
end

% The frist step will be to load all the corrmats into a common 3D array, then remove the missing nodes, then compute the louvain breakdown
%Create a matrix that will store all the:
%%% Assignments
Mall_reps = zeros(size(commonROIs,1), size(subList, 2));

%%% Start the subject forloop
tmp = 1;
for subject = subList;
    if strcmp(subject,'sub-029') == 1 || strcmp(subject,'sub-033') == 1;
        continue
    end
    
    if strcmp(atlas, 'Schaefer400node')
        corrmatFname = [subdir char(subject) '/' char(subject) '_task-rest_run-BOTH_space-MNI152NLin2009cAsym_desc-residuals_variant-24p_Acc6_009_1_Schaefer400node_mean_corrmat.npy'];
    elseif strcmp(atlas, 'brainnetome')
        corrmatFname = [subdir char(subject) '/' char(subject) '_task-rest_run-BOTH_space-MNI152NLin2009cAsym_desc-residuals_variant-24p_Acc6_009_1_brainnetome_mean_corrmat.npy'];
    end
    
    corrmat = readNPY(corrmatFname);
    % Load the subjects correlation matrix
    % Remove all the missing ROIs rows and columns
    
    for mN = flipud(missingROIs+1) % Make sure to add 1 because this list comes from Python with its 0-indexing
        corrmat(mN,:) = [];
        corrmat(:,mN) = [];
    end
    
    % Zero out all the diagonal entries since BCT states that there should be no self-connections
    corrmat_noZD = corrmat - diag(diag(corrmat));
    
    % Zero out all negative edges entirely
    corrmat_noZD(find(corrmat_noZD<0))=0;
    
    
    % 1. Compute modularity, including all positive weights
    %%%%%%%% Maybe include an option to play with the gamma parameter %%%%%%%%%%%%%%%
    for i = 1:Nmods
        [Mall, Qall] = community_louvain(corrmat_noZD, '', '', 'modularity');
        Mall_reps(:,i) = Mall;
    end
    
    % 2. Create an agreement matrix (a proportion, (0-N)/N, of times a node is in the same module)
    Dall = agreement(Mall_reps)/Nmods;
    
    % 3. Get consensus partition, running Nreps times, with a threshold of tau
    CIU_all = consensus_und(Dall, tau, Nreps);
    
    % 4. Calculate modularity given the partition
    [Mall, Qall] = community_louvain(corrmat_noZD, '', CIU_all, 'modularity');
    
    
    Qall_all(tmp, 1) = str2num(subject{:}(5:7));
    Qall_all(tmp, 2) = [Qall];
    
    tmp = tmp + 1;
    
    fileStr = ['Louvain_ModVals_FD50excld_ConsenPart_allPos_' atlas];
    
    save([outdir fileStr '.mat'], 'Qall_all')
    csvwrite([outdir fileStr '.csv'], Qall_all)
end