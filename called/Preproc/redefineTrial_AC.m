
function data = redefineTrial_AC(data, TOI)

% Define time and Time indices
Time = data.time{1};
TInd = find(Time >= TOI(1) & Time <= TOI(2));

% Loop over trials and adjust vectors

NTrials = length(data.trial);
Trials = data.trial;
data = rmfield(data, 'trial');

for tt = 1:NTrials
    data.trial{tt} = Trials{tt}(:,TInd);
    data.time{tt} = data.time{tt}(TInd);
end

data = rmfield(data, 'sampleinfo');