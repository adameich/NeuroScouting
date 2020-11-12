% 1. Compute modularity N times and save all output in an NNodes x N mod reps matrix
% 2. Calculate modularity given the partition
% 3. Save output partitions, modularity, and number modules

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

% Load in the Yeo 7-Network partition affiliations
if strcmp(atlas, 'Schaefer400node')
    YeoAff = load('/Users/Eichenbaum/HWNI/Experiments/nScout/data/resources/parcellations/Schaefer2018_LocalGlobal/LG400_7Network_memberships.txt');
    missingROIsFName = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/analyses/missingROIs_list_thresh-0.25_Schaefer400node.npy';
    commonROIsFName = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/analyses/commonROIs_list_thresh-0.25_Schaefer400node.npy';
    missingROIs = readNPY(missingROIsFName);
    commonROIs = readNPY(commonROIsFName);
    
elseif strcmp(atlas, 'brainnetome')
    YeoAff = load('/Users/Eichenbaum/HWNI/Experiments/nScout/data/resources/parcellations/brainnetome/brainnetome_7Network_memberships.txt');
    missingROIsFName = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/analyses/missingROIs_list_thresh-0.25_Brainnetome.npy';
    commonROIsFName = '/Users/Eichenbaum/HWNI/Experiments/nScout/data/analyses/commonROIs_list_thresh-0.25_Brainnetome.npy';
    missingROIs = readNPY(missingROIsFName);
    commonROIs = readNPY(commonROIsFName);
    %%% Brainnetome has subcortical nodes, Yeo does not. Make sure to add
    %%% the subcortical nodes as "missing" if using a predefinied partition
    subcortNodes = [211:246] - 1; % Going to add 1 to them later because the MissingNode numpy file is python 0-indexed, while the Yeo .txt file is matlab 1-indexed
end

for mN = flipud(missingROIs+1) % Make sure to add 1 because this list comes from Python with its 0-indexing
    YeoAff(mN) = []; % Remove the missing nodes from the Yeo partition file's index
end

subList = get_sublist(subdir);

if testBug;
    subList = {'sub-201'};
end

% The frist step will be to load all the corrmats into a common 3D array, then remove the missing nodes, then compute the louvain breakdown
%Create a matrix that will store all the:
%%% Assignments

%%% Start the subject forloop
tmp = 1;
for subject = subList;
    if strcmp(subject,'sub-029') == 1 || strcmp(subject,'sub-033') == 1; % Remove high in-scanner motion subjects
        continue
    end
    
    if strcmp(atlas, 'brainnetome')
        % Load the subjects correlation matrix
        corrmatFname = [subdir char(subject) '/' char(subject) '_task-rest_run-BOTH_space-MNI152NLin2009cAsym_desc-residuals_variant-24p_Acc6_009_1_brainnetome_mean_corrmat.npy'];
        corrmat = readNPY(corrmatFname);
        
        %%% nodes to remove from corrmat
        nodesToRemove = missingROIs;
        nodesToRemove(length(nodesToRemove)+1 : length(nodesToRemove) + length(subcortNodes)) = subcortNodes;
        % Remove all the missing ROIs rows and columns
        for mN = fliplr(nodesToRemove+1) % Make sure to add 1 because this list comes from Python with its 0-indexing
            corrmat(mN,:) = [];
            corrmat(:,mN) = [];
        end
    elseif strcmp(atlas, 'Schaefer400node')
        % Load the subjects correlation matrix
        corrmatFname = [subdir char(subject) '/' char(subject) '_task-rest_run-BOTH_space-MNI152NLin2009cAsym_desc-residuals_variant-24p_Acc6_009_1_Schaefer400node_mean_corrmat.npy'];
        corrmat = readNPY(corrmatFname);
        
        % Remove all the missing ROIs rows and columns
        for mN = flipud(missingROIs+1) % Make sure to add 1 because this list comes from Python with its 0-indexing
            corrmat(mN,:) = [];
            corrmat(:,mN) = [];
        end
    end
    
    %%% Just for fun, stack all the mean corrmats
    corrmat_all(:,:, tmp) = corrmat;
    
    % Zero out all the diagonal entries since BCT states that there should be no self-connections
    corrmat_noZD = corrmat - diag(diag(corrmat));
    
    % Zero out all negative edges entirely
    corrmat_noZD(find(corrmat_noZD<0))=0;
    
    [Q1] = compute_modularity(corrmat_noZD, [], YeoAff, 'modularity');
    
    Qall_all(tmp, 1) = str2num(subject{:}(5:7));
    Qall_all(tmp, 2) = [Q1];
    
    fileStr = ['Louvain_ModVals_FD50excld_Yeo7NetPart_allPos_' atlas];
    
    save([outdir fileStr '.mat'], 'Qall_all')
    csvwrite([outdir fileStr '.csv'], Qall_all)
    
    %%% Compute Global Efficiency %%%
    Eglob = efficiency_bin(corrmat_noZD);
    
    GE(tmp, 1) = str2num(subject{:}(5:7));
    GE(tmp, 2) = [Eglob];
    fileStrGE = ['Louvain_GEVals_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrGE '.mat'], 'GE')
    csvwrite([outdir fileStrGE '.csv'], GE)
    
    %%% Compute Participation Coefficient %%%
    P = participation_coef(corrmat_noZD, YeoAff);
    
    %%% Mean Brain-wide PC
    mean_pc = mean(P);
    
    %%% Mean PC of Visual Network
    Vis_pc = mean(P(find(YeoAff==1)));
    
    %%% Mean PC of Somatosensory-Motor Network
    Mot_pc = mean(P(find(YeoAff==2)));
    
    %%% Mean PC of Dorsal Attention Network
    DAN_pc = mean(P(find(YeoAff==3)));
    
    %%% Mean PC of Ventral Attention / Salience Network
    SalVAN_pc = mean(P(find(YeoAff==4)));
    
    %%% Mean PC of Limbic Network
    Limb_pc = mean(P(find(YeoAff==5)));
    
    %%% Mean PC of Control Network
    Cont_pc = mean(P(find(YeoAff==6)));
    
    %%% Mean PC of Deafult Mode Network
    DMN_pc = mean(P(find(YeoAff==7)));
    
    %%% Mean PC of Everything Else Network
    allElse_pc = mean(P(find(YeoAff>2)));
    
    PallwholeBrainMean(tmp, 1) = str2num(subject{:}(5:7));
    PallwholeBrainMean(tmp, 2) = [mean_pc];
    fileStrwholeBrainMean = ['Louvain_PCVals-wholeBrainMean_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrwholeBrainMean '.mat'], 'PallwholeBrainMean')
    csvwrite([outdir fileStrwholeBrainMean '.csv'], PallwholeBrainMean)
    
    PallVis(tmp, 1) = str2num(subject{:}(5:7));
    PallVis(tmp, 2) = [Vis_pc];
    fileStrVis = ['Louvain_PCVals-Vis_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrVis '.mat'], 'PallVis')
    csvwrite([outdir fileStrVis '.csv'], PallVis)
    
    PallMot(tmp, 1) = str2num(subject{:}(5:7));
    PallMot(tmp, 2) = [Mot_pc];
    fileStrMot = ['Louvain_PCVals-Mot_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrMot '.mat'], 'PallMot')
    csvwrite([outdir fileStrMot '.csv'], PallMot)
    
    PallDAN(tmp, 1) = str2num(subject{:}(5:7));
    PallDAN(tmp, 2) = [DAN_pc];
    fileStrDAN = ['Louvain_PCVals-DAN_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrDAN '.mat'], 'PallDAN')
    csvwrite([outdir fileStrDAN '.csv'], PallDAN)
    
    PallSalVAN(tmp, 1) = str2num(subject{:}(5:7));
    PallSalVAN(tmp, 2) = [SalVAN_pc];
    fileStrSalVAN = ['Louvain_PCVals-VANSal_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrSalVAN '.mat'], 'PallSalVAN')
    csvwrite([outdir fileStrSalVAN '.csv'], PallSalVAN)
    
    PallLimb(tmp, 1) = str2num(subject{:}(5:7));
    PallLimb(tmp, 2) = [Limb_pc];
    fileStrLimb = ['Louvain_PCVals-Limb_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrLimb '.mat'], 'PallLimb')
    csvwrite([outdir fileStrLimb '.csv'], PallLimb)
    
    PallCont(tmp, 1) = str2num(subject{:}(5:7));
    PallCont(tmp, 2) = [Cont_pc];
    fileStrCont = ['Louvain_PCVals-Cont_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrCont '.mat'], 'PallCont')
    csvwrite([outdir fileStrCont '.csv'], PallCont)
    
    PallDMN(tmp, 1) = str2num(subject{:}(5:7));
    PallDMN(tmp, 2) = [DMN_pc];
    fileStrDMN = ['Louvain_PCVals-DMN_FD50excld_Yeo7NetPart_allPos_' atlas];
    %         save([outdir fileStrDMN '.mat'], 'PallDMN')
    csvwrite([outdir fileStrDMN '.csv'], PallDMN)
    
    PallElse(tmp, 1) = str2num(subject{:}(5:7));
    PallElse(tmp, 2) = [allElse_pc];
    fileStrElse = ['Louvain_PCVals-Else_FD50excld_Yeo7NetPart_allPos_' atlas];
    save([outdir fileStrElse '.mat'], 'PallElse')
    csvwrite([outdir fileStrElse '.csv'], PallElse)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% [Visual Network] %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Within-Network Mean Connectivity %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VisNetW = corrmat_noZD(find(YeoAff==1), find(YeoAff==1));
    VisNetW_all(:,:,tmp) = VisNetW;
    VisNetWvec = VisNetW(find(~tril(ones(size(VisNetW)))));

    meanVisNetW = mean(VisNetWvec);
    
    meanVisW_all(tmp, 1) = str2num(subject{:}(5:7));
    meanVisW_all(tmp, 2) = [meanVisNetW];
    fileStrVisW = ['Louvain_meanFCVals-VisW_FD50excld_Yeo7NetPart_allPos_' atlas];
    csvwrite([outdir fileStrVisW '.csv'], meanVisW_all)
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% [Motor Network] %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Within-Network Mean Connectivity %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MotNetW = corrmat_noZD(find(YeoAff==2), find(YeoAff==2));
    MotNetW_all(:,:,tmp) = MotNetW;
    MotNetWvec = MotNetW(find(~tril(ones(size(MotNetW)))));
    
    meanMotNetW = mean(MotNetWvec);
    
    meanMotW_all(tmp, 1) = str2num(subject{:}(5:7));
    meanMotW_all(tmp, 2) = [meanMotNetW];
    fileStrMotW = ['Louvain_meanFCVals-MotW_FD50excld_Yeo7NetPart_allPos_' atlas];
    csvwrite([outdir fileStrMotW '.csv'], meanMotW_all)
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% [Dorsal Attention Network] %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Within-Network Mean Connectivity %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DANNetW = corrmat_noZD(find(YeoAff==3), find(YeoAff==3));
    DANNetW_all(:,:,tmp) = DANNetW;
    DANNetWvec = DANNetW(find(~tril(ones(size(DANNetW)))));
    
    meanDANNetW = mean(DANNetWvec);
    
    meanDANW_all(tmp, 1) = str2num(subject{:}(5:7));
    meanDANW_all(tmp, 2) = [meanDANNetW];
    fileStrDANW = ['Louvain_meanFCVals-DANW_FD50excld_Yeo7NetPart_allPos_' atlas];
    csvwrite([outdir fileStrDANW '.csv'], meanDANW_all)
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% [Ventral Attention / Salience Network] %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Within-Network Mean Connectivity %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VANSalNetW = corrmat_noZD(find(YeoAff==4), find(YeoAff==4));
    VANSalNetW_all(:,:,tmp) = VANSalNetW;
    VANSalNetWvec = VANSalNetW(find(~tril(ones(size(VANSalNetW)))));
    
    meanVANSalNetW = mean(VANSalNetWvec);
    
    meanVANSalW_all(tmp, 1) = str2num(subject{:}(5:7));
    meanVANSalW_all(tmp, 2) = [meanVANSalNetW];
    fileStrVANSalW = ['Louvain_meanFCVals-VANSalW_FD50excld_Yeo7NetPart_allPos_' atlas];
    csvwrite([outdir fileStrVANSalW '.csv'], meanVANSalW_all)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% [Limbic Network] %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Within-Network Mean Connectivity %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LimbNetW = corrmat_noZD(find(YeoAff==5), find(YeoAff==5));
    LimbNetW_all(:,:,tmp) = LimbNetW;
    LimbNetWvec = LimbNetW(find(~tril(ones(size(LimbNetW)))));
    
    meanLimbNetW = mean(LimbNetWvec);
    
    meanLimbW_all(tmp, 1) = str2num(subject{:}(5:7));
    meanLimbW_all(tmp, 2) = [meanLimbNetW];
    fileStrLimbW = ['Louvain_meanFCVals-LimbW_FD50excld_Yeo7NetPart_allPos_' atlas];
    csvwrite([outdir fileStrLimbW '.csv'], meanLimbW_all)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% [Control Network] %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Within-Network Mean Connectivity %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ContNetW = corrmat_noZD(find(YeoAff==6), find(YeoAff==6));
    ContNetW_all(:,:,tmp) = ContNetW;
    ContNetWvec = ContNetW(find(~tril(ones(size(ContNetW)))));
    
    meanContNetW = mean(ContNetWvec);
    
    meanContW_all(tmp, 1) = str2num(subject{:}(5:7));
    meanContW_all(tmp, 2) = [meanContNetW];
    fileStrContW = ['Louvain_meanFCVals-ContW_FD50excld_Yeo7NetPart_allPos_' atlas];
    csvwrite([outdir fileStrContW '.csv'], meanContW_all)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% [Default Network] %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Within-Network Mean Connectivity %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DMNNetW = corrmat_noZD(find(YeoAff==7), find(YeoAff==7));
    DMNNetW_all(:,:,tmp) = DMNNetW;
    DMNNetWvec = DMNNetW(find(~tril(ones(size(DMNNetW)))));
    
    meanDMNNetW = mean(DMNNetWvec);
    
    meanDMNW_all(tmp, 1) = str2num(subject{:}(5:7));
    meanDMNW_all(tmp, 2) = [meanDMNNetW];
    fileStrDMNW = ['Louvain_meanFCVals-DMNW_FD50excld_Yeo7NetPart_allPos_' atlas];
    csvwrite([outdir fileStrDMNW '.csv'], meanDMNW_all)
    
    
    tmp = tmp + 1;
end