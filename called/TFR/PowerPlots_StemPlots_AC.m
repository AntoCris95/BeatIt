
function PowerPlots_AC(sess, BSLcorr)


% Set up

if nargin < 1
    sess = input('Session number?');
    BSLcorr= input('Use BSL-corrected data? 0= No, 1= Yes');
end

datatype= 0; disp('Using simple (0) data ');
datadir = getSessDir(sess);
datadir.output = fullfile(datadir.sess, 'DATAproc');

if BSLcorr == 1
    outdir= fullfile(datadir.output, 'TFR', 'BSLcorr', 'Results');
else
    outdir= fullfile(datadir.output, 'TFR', 'meancorr', 'Results');
end

load('C:\Users\p70068941\OneDrive\Work\UM\BAND\Projects\GeneralFunctions\TFR\mycolormap.mat')

% Param of interest
Labels= {'STD', 'STDs', 'STDw', 'DEV', 'DEVs', 'DEVw'};
Comparisons = {'STD s-w', 'DEV s-w', 'STD-DEV'};
PlannedContrasts(1,:) = [2,3];
PlannedContrasts(2,:) = [5,6];
PlannedContrasts(3,:) = [1,4];

FqOI = [1,40];
TimeOrig = [-1 1];
Amp(1,:) = [-1 1];
Amp(2,:) = [-3 3];
Amp(3,:) = [-2 2];


%% Here we go

for cond = 1:length(Comparisons)
    
    % Import data
    
    load(fullfile(outdir, ['GGA ' Labels{PlannedContrasts(cond,1)} '.mat']));
    data1 = data;
    
    load(fullfile(outdir, ['GGA ' Labels{PlannedContrasts(cond,2)} '.mat']));
    data2 = data; clear data
    
    % Define Time and freqs
    
    Time= linspace(TimeOrig(1),TimeOrig(2),size(data1,6));
    T_OI = [-.3 .3];
    Tind = find(Time >= T_OI(1) & Time <= T_OI(2));
    
    fqs = linspace(FqOI(1),FqOI(2),size(data1,5));
    [val, indi] = sort(fqs,'descend');
    
    % Power plots
    
    data1 = squeeze(data1);
    data2 = squeeze(data2);
    
    figure('Color', 'White');
    subplot(1,2,1)
    imagesc(squeeze(data1(end,indi,Tind)), Amp(cond,:)) %, [-1 1]
    colormap(mycolormap)
    colorbar
    
    a = gca; 
    a.FontSize = 14;
    a.FontWeight = 'bold';
    a.XTickLabel = round(Time(Tind(a.XAxis.TickValues)),1);
    a.YTickLabel = round(fqs(indi(a.YAxis.TickValues)));
    
    a.XLabel.String = 'Time (s)';
    a.XLabel.Color = [0 0 0];
    a.YLabel.String = 'Frequency (Hz)';
    a.YLabel.Color = [0 0 0];
    box off
    grid off
    
    title(Labels{PlannedContrasts(cond,1)})
    
    subplot(1,2,2)
    imagesc(squeeze(data2(end,indi,Tind)), Amp(cond,:)) %, [-1 1]
    colormap(mycolormap)
    colorbar
    
    a = gca; 
    a.FontSize = 14;
    a.FontWeight = 'bold';
    a.XTickLabel = round(Time(Tind(a.XAxis.TickValues)),1);
    a.YTickLabel = round(fqs(indi(a.YAxis.TickValues)));
    
    a.XLabel.String = 'Time (s)';
    a.XLabel.Color = [0 0 0];
    a.YLabel.String = 'Frequency (Hz)';
    a.YLabel.Color = [0 0 0];
    box off
    grid off
    
    title(Labels{PlannedContrasts(cond,2)})
    
    % Stem plots
    
    Colors = {'-b', '-r'};
    
    TOI1 = [-.2 -.05]; % pre-stim
    TOI2 = [0 .15]; % post-stim
    
    TimeInd1 = find(Time >= TOI1(1) & Time <= TOI1(2)); % pre-stim
    TimeInd2 = find(Time >= TOI2(1) & Time <= TOI2(2)); % post-stim
    
    datasel1 = [];
    datasel2 = [];
    
    % Delta
    
    FQ_Sel = 1:4; % delta
    FQ_Ind = find(fqs >= FQ_Sel(1) & fqs <= FQ_Sel(2));
    
    datasel1(1) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd1)));
    datasel2(1) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd1)));
    
    datasel1(3) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd2)));
    datasel2(3) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd2)));
    
    % Theta
    
    FQ_Sel = 4:8; % theta
    FQ_Ind = find(fqs >= FQ_Sel(1) & fqs <= FQ_Sel(2));
    
    datasel1(5) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd1)));
    datasel2(5) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd1)));
    
    datasel1(7) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd2)));
    datasel2(7) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd2)));
    
    % Alpha
    
    FQ_Sel = 8:12; % alpha
    FQ_Ind = find(fqs >= FQ_Sel(1) & fqs <= FQ_Sel(2));
    
    datasel1(9) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd1)));
    datasel2(9) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd1)));
    
    datasel1(11) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd2)));
    datasel2(11) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd2)));
    
    % Lowbeta
    
    FQ_Sel = 12:20; % lowbeta
    FQ_Ind = find(fqs >= FQ_Sel(1) & fqs <= FQ_Sel(2));
    
    datasel1(13) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd1)));
    datasel2(13) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd1)));
    
    datasel1(15) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd2)));
    datasel2(15) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd2)));
    
    % Highbeta
    
    FQ_Sel = 20:25; % lowbeta
    FQ_Ind = find(fqs >= FQ_Sel(1) & fqs <= FQ_Sel(2));
    
    datasel1(17) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd1)));
    datasel2(17) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd1)));
    
    datasel1(19) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd2)));
    datasel2(19) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd2)));
    
    % Gamma
    
    FQ_Sel = 25:40; % lowbeta
    FQ_Ind = find(fqs >= FQ_Sel(1) & fqs <= FQ_Sel(2));
    
    datasel1(21) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd1)));
    datasel2(21) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd1)));
    
    datasel1(23) = nanmean(nanmean(data1(end,FQ_Ind,TimeInd2)));
    datasel2(23) = nanmean(nanmean(data2(end,FQ_Ind,TimeInd2)));
    
    datasel1 = [0 0 datasel1 0 0];
    datasel2 = [0 0 datasel2 0 0];
    
    % Plotting
    
    figure('color', 'white');
    stem(datasel1(1:2:end), 'filled', Colors{1}, 'linewidth', 2); hold on
    stem((1:round(length(datasel2)/2))+0.2,datasel2(1:2:end), 'filled', Colors{2}, 'linewidth', 2);
    
    % Take out 0 points
    % Beginning
    stem(datasel1(1), 'filled', '-w', 'linewidth', 2);
    stem(1+0.2,datasel2(1), 'filled', '-w', 'linewidth', 2);
    % End
    stem(14, datasel1(end), 'filled', '-w', 'linewidth', 2);
    stem(14+0.2,datasel2(1), 'filled', '-w', 'linewidth', 2);
    
    % Last touch
    box off
    a = gca; a.XColor = [1 1 1];
    a.FontSize = 14;
    a.FontWeight = 'bold';
    a.YLabel.String = 'Norm %change';
    hline = refline([0,0]);
    hline.Color = 'k';
    h = legend;
    h.String{1} = Labels{PlannedContrasts(cond,1)};
    h.String{2} = Labels{PlannedContrasts(cond,2)};
    h.String(3:end) = [];
    legend boxoff
    
end

disp('Remember to add amplitude bars')



