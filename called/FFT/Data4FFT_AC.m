
%% Pipeline for Fourier analyses

% This function uses partially preprocessed data from Preprocess_AC
% then creates two Data structures for the isochronous and random conditions
% containing multiple segments of 20s (Segmentation4FFT_AC).
% Next, it undergoes discrete fourier transform by invoking goFourier_AC,
% creates GA per conditions and plots the data with shaded error bars.


function Data4FFT_AC(sess)


if nargin < 1
    sess = input('Session number?');
end

% Prepare data

datadir = getSessDir(sess);
outdir = fullfile(datadir.sess, 'DATAproc'); %, 'TFR');
outdir2 = fullfile(outdir, 'FFT');

if ~exist(outdir2)
    mkdir(outdir2)
end

% Set up

FQS = [0.2 30]; % for fft spectrum
Tones_OI = [3 13];

load(fullfile(datadir.main, 'chLabels.mat'));
chan_labs = Labels; clear Labels

Ch2Rem = {'EOGV', 'EOGH', 'A1', 'A2'};
indch = find(ismember(chan_labs,Ch2Rem) < 1);
chan_labs = chan_labs(indch);

% Channels of interest
chLabels = 'E:\DATA\Humans - Beatit\DATAproc\chs_OI.mat';
load(chLabels); %chs_OI

indch = [];
for cc = 1:length(chs_OI)
    inds = find(strcmp(chan_labs, chs_OI{cc}));
    indch = [indch inds];
end


%% Here we go

files = dir(fullfile(outdir, 'ICA *.mat'));

for ff = 1:length(files)
    
    fprintf('\n Looping over file %d / %d', ff, length(files))
    
    saveName1 = ['ep4FFT ' files(ff).name];
    saveName2 = ['FFT ' saveName1];
    saveName1 = fullfile(outdir2, saveName1);
    saveName2 = fullfile(outdir2, saveName2);
    
    if ~exist(saveName2)
        if ~exist(saveName1)
            
            load(fullfile(outdir, files(ff).name));
            
            fprintf('\n Setting up Segmentation for FFT')
            data = Segmentation4FFT_AC(data, datadir, ff);
            save(saveName1, 'data');
            
        else
            load(saveName1)
        end
        
        % Prepare for fft analyses
                
        Time = data.time{1};
        T_OI = [Tones_OI(1)*data.ISO-data.ISO Tones_OI(2)*data.ISO-data.ISO];
        TInd = find(round(Time,3)>= T_OI(1) & round(Time,3)<= T_OI(2));
        
        fprintf('\n Performing FFT')
        data = goFourier_AC(data, FQS, indch, TInd);
        %         data = goFourier2_AC(data, FQS, indch, TInd); % avg trials
        save(saveName2, 'data');
        
    end
end

clear data files

% Create GAs

fprintf('\n Preparing GAs')

files1 = dir(fullfile(outdir2, 'FFT *.mat'));

% Load in FFT data

ind = 2; % from the second fq
Data1 = []; % normal
Data2 = []; % mean correction with neighbor freq

for ff = 1:length(files1)
    load(fullfile(outdir2, files1(ff).name));
    Data1(ff,:) = nanmean(data.fft(:,ind:end),1);
    Data2(ff,:) = nanmean(data.fft2(:,ind:end),1);
end


% Compute avg

ISO= [];
ISO.xaxis= data.FqOI(1:size(Data1,2));
% ISO.xaxis= data.FqOI;
ISO.nanmean= nanmean(Data1);
ISO.SE= (std(Data1, 'omitnan'))/sqrt(size(Data1,1));
ISO.nanmean2= nanmean(Data2);
ISO.SE2= (std(Data2, 'omitnan'))/sqrt(size(Data2,1));

saveName = fullfile(outdir2, 'GA FFT.mat');

save(saveName, 'ISO');


%% Plot

fprintf('\n Ready for Single-subj Plotting')

figure('Color', 'White'); hold on, box off
n = 1;

for ff = 1:2:size(Data1, 1)
    try
        subplot(4,5,n); hold on
        plot(ISO.xaxis, Data1(ff:ff+1,:));
        title(sprintf('FFT Subj%d', n))
        n = n+1;
    end
end

figure('Color', 'White'); hold on, box off
n = 1;

for ff = 1:2:size(Data2, 1)
    try
        subplot(4,5,n); hold on
        plot(ISO.xaxis, Data2(ff:ff+1,:));
        title(sprintf('FFT Subj%d fqcorr', n))
        n = n+1;
    end
end

fprintf('\n Ready for GA Plotting')

% Normal

figure('Color', 'White'); hold on, box off
shadedErrorBar(ISO.xaxis, ISO.nanmean, ISO.SE, 'lineProps','-b'); %'lineProps',
title(['GA, FFT ISO'])
a = gca;
a.XAxis.TickValues = 0 : 1 : 30;

% Fq corr
figure('Color', 'White'); hold on, box off
shadedErrorBar(ISO.xaxis, ISO.nanmean2, ISO.SE2, 'lineProps','-b'); %'lineProps',
title(['GA, FFT ISO fqcorr'])
a = gca;
a.XAxis.TickValues = 0 : 1 : 30;
a.FontSize = 14;
a.FontWeight = 'bold';

xlabel('\fontsize{14} Frequency (Hz)', 'fontweight', 'bold')
ylabel('\fontsize{14} Normalized Pow', 'fontweight', 'bold')

