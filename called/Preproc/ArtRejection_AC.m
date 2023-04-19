

%% Artifact rejection based on amplitude value

% This function performs artifact rejection by means of FieldTrip toolbox
% Input data is the output of ft_preprocessing.



function [data, Nartifacts] = ArtRejection_AC(data)

if length(size(data.trial{1})) > 2
    error('ArtRejection only works on 2D data')
end

NTrials = length(data.trial);
try Events = data.events; end

% Set up thresholds

MaxVal = 100; % before was 85uV
MinVal = -MaxVal;

% Select windows of interest for artifact detection

Time = data.time{1};
TWind = [-.1 .3];
TIndi = find(Time >= TWind(1) & Time <= TWind(2)); 

artifact = [];
for tt = 1:length(data.trial)
    for cc = 1:size(data.trial{tt},1)
        if find(abs(data.trial{tt}(cc,TIndi)) >= MaxVal)
            artifact(end+1) = tt;
        end
    end
end

artifact = unique(artifact);
Nartifacts = length(artifact);

% Remove artifacts

if ~isempty(artifact)
    data.sampleinfo(artifact,:) = [];
    data.trial(artifact) = [];
    data.time(artifact) = [];
    try
    data.events(artifact) = [];
    end
end
