

% % Artifact suppresion (AS) by mean of temporal linear interpolation % %


% We will perform AS based on channel-by-channel and trial-by-trial
% estimations of outliers by inspection of mean/median and SD values.
% We work on a temporary data structure and return artifact-clean data to
% the orig data to continue preprocessing.

% We make use of the parallel computing toolbox to speed up computations


% Created in July 2020
% Edited in October 2020
% Written by Antonio Criscuolo


function [data, Artifacts] = ArtSuppression_AC(data, Normalize, TOI)

   
% Make sure there are no complex exponentials

if ~isreal(data.trial{1})
    disp('ICA did not work out. Interrupting')
    return
end

% Get back to a data matrix %

Ntrials= length(data.trial);
Nchannels= length(data.label);
Ntpoints= length(data.time{1});
try Events = data.events; 
SampleInfo = data.sampleinfo; end

allData= zeros(Ntrials, Nchannels, Ntpoints); %pre-allocate

for nt = 1:Ntrials
    allData(nt,:,:)= data.trial{nt};
end

% Normalize data 
if Normalize
    [normData] = normalizeData_AC(allData,1);
    Noise = [3 4];
else
    normData = allData; % when running AS with normal amplitude values
    Noise = [80 100];
end

% Time of interest
if ~exist('TOI', 'var')
    TOI = [-.4 .4];
%     TOI = [data.time{1}(1) data.time{1}(end)];
end
TInd = find(data.time{1} >= TOI(1) & data.time{1} <= TOI(2));

% % Initiate parallel loops % %

% Iterate through channels %
% Estimate a threshold value based on median/mean and SD %
% and perform artifact suppression %

parfor cc = 1:Nchannels
   [allData(:,cc,:), Artifacts{cc}] = newParallelLoops_ArtSuppression_AC(data, allData(:,cc,:), cc, normData(:,cc,:), TInd, Noise);
end

% Get back to a data structure for FieldTrip

for nt = 1:size(allData,1)
    data.trial{nt}= squeeze(allData(nt,:,:));
end

try data.events = Events; 
data.sampleinfo = SampleInfo; end

