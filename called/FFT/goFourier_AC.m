
function data = goFourier_AC(data, FQS, Chs, TInd)

% This function handles 3D data structures with Fieldtrip format 
% and performs 2D Fourier decomposition looping over trials and channels. 
% Next, data is resampled, normalized and amplitude-corrected at 
% every frequency-bin with the mean of neighboring frequency-bins. 

% Use as follows:
% data = goFourier_AC(data, FQS, Chs)

% with data formatted in Fieldtrip format;
% FQS as 1x2 vector of frequencies of interest (e.g. [1 30]);
% Chs as 'all' in case you want to use all channels or
% as a 1xN vector with channel indices. 

% At the bottom, there are some plotting functions to edit if needed. 

% Created in April 2021
% Written by Antonio Criscuolo


% Set up

if nargin < 2
    FQS = [0.2 30];
    Chs = 'all';
    TInd = 1:length(data.trial{1});
elseif nargin >= 2
    if ~exist('Chs', 'var')
        Chs = 'all';
    end
    if ~exist('FQS', 'var')
        FQS = [0.2 30];
    end
    if ~exist('TInd', 'var')
        TInd = 1:length(data.trial{1});
    end
end

if strcmp(Chs, 'all')
    Chs = 1:length(data.label);
else
    try % specific labels 
        Chs = find(ismember(data.label, Chs));
    end % or chs N?
end
if isfield(data, 'chSuppr')
    indi = ismember(Chs, unique([data.chSuppr{1}', data.chSuppr{2}']));
    Chs(indi) = [];
end

% Make sure there are no nans
% Interpolate nans in time
% Focus on chs of interest

SampleInfo = data.sampleinfo;
data = removefields(data, 'sampleinfo');

cfg = [];
cfg.channel = data.label(Chs);
data = ft_selectdata(cfg, data);

for tr = 1:length(data.trial)
    data.trial{tr} = fillmissing(data.trial{tr},'pchip',2);
end

UseFT = 1;
Fs = data.fsample; % sampling Fq
L = length(data.time{1}(TInd))-1; % length 

% Use fourier transform

if UseFT 
    
    % Use FTrip implementation of fft
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.output = 'fourier';
    cfg.pad = 'nextpow2';
    cfg.foi = FQS(1):.01:FQS(2);
    % cfg.keeptrials = 'yes';
    fr = ft_freqanalysis(cfg, data);
    fr.abs_fourierspctrm = abs(fr.fourierspctrm);
    
    Data = squeeze(nanmean(fr.abs_fourierspctrm));
    
    % 1. Power is the abs of the fourier spectrum
    % Amplitude is the sqrt of the fft ouput. 
    % Note that FTrip is already giving you power. 
    Y = Data; 
    
else
    
    % Use matlab-based fft
    for tt = 1:length(data.trial)
        tempdata(tt,:) = fft(nanmean(data.trial{tt}(Chs,TInd)), [],2);
    end
    
    % 1. Power is the abs of the fourier spectrum, 
    % which is expressed in magnitude. 
    % Amplitude (in volts) is the sqrt of the fft ouput. 
    tempdata = abs(tempdata);
    
    % Correct the length and amplitude
    Data = tempdata(:,1:round(L/2))*2;
    
    % Normalize
%     Data = Data ./norm(Data, 'fro');
    Data = Data ./ L;
%     Data = nanmean(Data);
    
    % Select fq of interest
    Y = Data;
    
end


% Data = sqrt(Data); % move from power to amplitude


% 2. Resample (if required)

ResampFact = 1;

if ResampFact > 1
    Y1 = resample(Y,1,ResampFact); % resample at specified factor 
    
    % Adjust parameters after resampling
    L = length(Y1)-1;
    Fs = Fs/ResampFact;
    frequencies = 0:1/(L/Fs):Fs/2; 
else
    Y1 = Y;
end


% 3. Normalize

% maybe not necessary?
% Y2 = normalizeData_AC(Y1);
Y2 = Y1;

% 3. Correct for neighboring fqs

for tt = 1:length(Y2)
    
    prepad = 4; % 4 neighboring fq
    postpad = 4; 
    
    % shorten if at the edges
    if prepad >= tt 
        prepad = tt-1;
    elseif (length(Y2) - tt) <= postpad
        postpad = length(Y2) - tt;
    end
    
    Indi = tt-prepad : tt+postpad;
    
    meanNeigh = nanmean(Y2(:,Indi),2);
    Y3(:,tt) = abs(Y2(:,tt) - meanNeigh);
    
end


% Send back data

try
    frequencies = fr.freq;    
catch
    frequencies= [0:1/(L/Fs):Fs/2];
end

FqOI = FQS;
fqInd = find(frequencies >= FqOI(1) & frequencies <= FqOI(end));

clear data
data.name = 'fft resampled normalized meancorr';
data.fsample = Fs/ResampFact;
data.FqOI = frequencies(fqInd);
data.fft = single(Y2(:,fqInd)); 
data.fft2 = single(Y3(:,fqInd));


% Plot

% figure('Color', 'White');
% subplot(4,1,1)
% plot(frequencies(1:300),Y(1:300))
% title('fft')
% 
% subplot(4,1,2)
% plot(NewFq(fqInd),Y1(:,fqInd))
% title('+resampled')
% 
% subplot(4,1,3)
% plot(NewFq(fqInd),Y2)
% title('+normalized')
% 
% subplot(4,1,4)
% plot(NewFq(fqInd),Y3)
% title('+meancorr')
% 
% % Plot final
% 
% figure('Color', 'White');
% plot(NewFq(fqInd),Y3, 'LineWidth',2)
% box off
% title('Fourier spectrum')



