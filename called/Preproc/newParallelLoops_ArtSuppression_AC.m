

% Parallel loop implementation of
% % Artifact suppresion (AS) by mean of temporal cubic interpolation % %


% We will perform AS based on channel-by-channel and trial-by-trial
% estimations of outliers by inspection of mean/median and SD values.
% We work on a temporary data structure and return artifact-clean data to
% the orig data to continue preprocessing.

% We make use of the parallel computing toolbox to speed up computations


% Created in October 2020
% Edited in January 2021
% Written by Antonio Criscuolo


function [allData, Artifacts] = newParallelLoops_ArtSuppression_AC(data, allData, cc, normData, TInd, Noise)


% Set up

Ntrials = length(data.trial);
chan_labs = data.label;
MaxArtTime = numel(normData(1,1,TInd))/3; % Max length of artifact windows

% Identify max and min values

if size(normData,1) > 30 % if there is a decent amount of trials
    
    SDfactor= 2; % adjust according to data
    meanData1 = nanmean(abs(normData(abs(normData) > .1))); % mean of non-zero data values
    SDdata= nanmean(std(normData(abs(normData) > meanData1))); % SD of out-of-mean values
    MaxVal = round(meanData1 + SDfactor*SDdata,1);
    
    if ~(MaxVal <= Noise(2) &&  MaxVal >= Noise(1)) % if it is not within a specified range, set it to a standard value
        MaxVal = mean(Noise);
    end
    
else
    MaxVal =  mean(Noise);
end
MinVal = - MaxVal;

% Create temporary data structure to allow
% channel-by-channel computations
% NB: Fieldtrip doesn't seem to like single-ch operations, so
% will try and trick it

tempdata = [];
tempdata.fsample = data.fsample;
tempdata.label{1} = chan_labs{cc};
tempdata.label{2} = chan_labs{cc};

% Use normalize data in time-windows of interest

for ttt = 1:length(data.trial)
    tempdata.time{ttt} = data.time{1}(TInd);
    tempdata.trial{ttt}(1,:) = normData(ttt,1,TInd); % work on normalized data
    tempdata.trial{ttt}(2,:) = normData(ttt,1,TInd); % will replace with orig data later
end

% Go with identification of outliers %

artifact = [];
newartifact = [];

cfg = [];
cfg.artfctdef.threshold.channel = {'all'};  %; chan_labs{cc+2}
cfg.artfctdef.threshold.bpfilter = 'no';
cfg.artfctdef.threshold.min = MinVal; % only marks time-points exceeding thresh
cfg.artfctdef.threshold.max = MaxVal;

[cfg, artifact] = ft_artifact_threshold(cfg, tempdata);


% Get back to non-normalized data

tempdata = rmfield(tempdata,'trial');

for ttt = 1:length(data.trial)
    
    % non-norm data
    tempdata.trial{ttt}(1,:) = allData(ttt,1,TInd);
    tempdata.trial{ttt}(2,:) = allData(ttt,1,TInd);
    
end

% Add padding to artifact time-windows

if ~isempty(artifact)
    
    padding = 0.05; % padding - 2use on both pre- and post-windows
    tpoints= round(padding / (1/data.fsample)); % N of time-points for the shift
    
    artifact(:,1)= artifact(:,1) - tpoints; % add pre-
    artifact(:,2)= artifact(:,2) + tpoints; % add post-
    
    % Avoid negative nums
    if find(artifact <= 0)
        indi = find(artifact <= 0);
        artifact(indi) = artifact(indi) + tpoints+1;
    end
    
    % Calculate length of artifact time-points
    for trl = 1:size(artifact,1)
        ToT = artifact(trl,2) - artifact(trl,1);
        
        if ToT < MaxArtTime % check we have a min of time-points
            newartifact(end+1,:) = artifact(trl,:);
        end
    end
    
    artifact = newartifact; clear newartifact
end

% % Replace time-points with artificats with NaNs % %

if ~isempty(artifact)
    
    cfg = [];
    cfg.artfctdef.reject = 'nan'; % replace artifacts with nans
    cfg.artfctdef.minaccepttim = 0.1;
    cfg.artfctdef.visual.artifact = artifact;
    
    tempdata = ft_rejectartifact(cfg, tempdata);
    
end


% Get back to full time vector

tempData = [];
tempData.trial = tempdata.trial;
tempdata = rmfield(tempdata,'trial');

for ttt = 1:length(data.trial)
    
    tempdata.time{ttt} = data.time{1}(TInd);
    
    % non-norm data
    tempdata.trial{ttt}(1,:) = allData(ttt,1,TInd);
    tempdata.trial{ttt}(2,:) = allData(ttt,1,TInd);
    
    % replace TInd with AS data
    tempdata.trial{ttt}(1,TInd) = tempData.trial{ttt}(1,:);
    tempdata.trial{ttt}(2,TInd) = tempData.trial{ttt}(1,:);
    
end


% % Use temporal interpolation to replace NaNs and keep (almost) all trials % %

try % if there are overlapping time-windows with nans, it won't work
    cfg = [];
    cfg.prewindow = 0.5; % if possible, uses interpolation with 500ms before+after artefact, otherwise it shortens the window automatically
    cfg.postwindow = 0.5;
    cfg.method = 'pchip';
    % 'linear' - default, linear interpolation
    % 'pchip' - shape-preserving piecewise cubic interpolation
    % 'v5cubic' - cubic interpolation from Matlab 5
    tempdata = ft_interpolatenan(cfg, tempdata);
end

% Another round of interpolation in case of missing values

clear tempData
tempData = zeros([length(tempdata.trial) length(tempdata.time{1})]);
for trl = 1:length(tempdata.trial)
    tempData(trl,:) = tempdata.trial{trl}(1,:);
end

try % if there are overlapping time-windows with nans, it won't work
    tempData= fillmissing(tempData, 'pchip'); %     'spline', 'pchip', 'makima'
end

% % Back to the full data structure % %

for trl = 1:size(tempData,1)
    allData(trl,1,:) = tempData(trl,:);
end


% % Artifacts summary

Nartifacts = size(artifact, 1) / 2; % divided by 2 because we replicate the same ch. so, actually it's only half of the artifacts

Artifacts.Num= Nartifacts;
if isnan(Artifacts.Num)
    Artifacts.Num = 0;
end
Artifacts.ThreshVal= MaxVal;
Artifacts.Length = artifact(:,2) - artifact(:,1);

% Trace back time-windows with artifacts

TrialSample(:,1) = 1:length(TInd):length(TInd)*length(tempdata.trial);
TrialSample(:,2) = length(TInd):length(TInd):length(TInd)*length(tempdata.trial);

TrlArt = [];
for tt = 1:size(artifact,1)

    ind = find(artifact(tt,1) >= TrialSample(:,1) & artifact(tt,1) <= TrialSample(:,2));
    TrlArt = [TrlArt ind];

    tmpTimeWin = TrialSample(ind,1) : TrialSample(ind,2);
    TimeWin = find(ismember(tmpTimeWin, [artifact(tt,1):artifact(tt,2)]));
    TimeWin = tempdata.time{1}(TimeWin);

    Artifacts.TimeWindows{tt} = TimeWin;
end

TrlArt = unique(TrlArt);
Artifacts.Trials = TrlArt;

% clear artifact


