
%% Plot time-frequency power plots

% Set up folder
path2load = 'E:\DATA\Humans - Beatit\DATAproc\TFR\meancorr\Results\GGA STD.mat';
load(path2load)

% Set up param
FqOI = [1,40];
fqs = linspace(FqOI(1),FqOI(2),size(data,5));
[val, indi] = sort(fqs,'descend');

TimeOrig = [-1 1];
Time= linspace(TimeOrig(1),TimeOrig(2),size(data,6));
T_OI = [-.2 .3];
Tind = find(Time >= T_OI(1) & Time <= T_OI(2));

% color map
load('C:\Users\p70068941\OneDrive\Work\UM\BAND\Projects\GeneralFunctions\TFR\mycolormap.mat')
Amp(1,:) = [-1 1];

% Go with plot

data = squeeze(data(:,:,:,:,indi,Tind));

figure('Color', 'White');
imagesc(squeeze(data(end,:,:)), Amp(1,:)); %, [-1 1]
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

% title('STD')

% Multiplots

path2load = 'E:\DATA\Humans - Beatit\chLabels.mat';
load(path2load)
chan_labs = Labels; clear Labels

%Rremove unwanted chs
Ch2Rem = {'EOGV', 'EOGH', 'A1', 'A2'}; 
indch = find(ismember(chan_labs,Ch2Rem) < 1);
chan_labs = chan_labs(indch);
% chs of interest
path2load = 'E:\DATA\Humans - Beatit\DATAproc\chs_OI.mat';
load(path2load);

data1 = [];
data1.label = chan_labs; 
data1.dimord = 'rpt_chan_freq_time';
data1.freq = fqs;
data1.time = Time(Tind);
data1.powspctrm = zeros([1, size(data)]); 
% data1.powspctrm(1,:,:,:) = data;
data1.powspctrm(:,end,:,:) = [];


cfg = [];
cfg.layout = 'C:\Users\p70068941\Documents\MATLAB\fieldtrip\template\layout\biosemi64.lay';
cfg.parameter = 'powspctrm';
cfg.highlight = 'on';
cfg.highlightchannel = chs_OI;
cfg.ylim = [-1 1];
cfg.colormap = mycolormap;
cfg.interpolatenan = 'no';

ft_topoplotTFR(cfg, data1)


%% Plot Beta Peaks as stem plots and violin plots

% Load-in Peak mat

path2load = 'E:\DATA\Humans - Beatit\DATAproc\TFR\meancorr\Results\TFR Peaks_MainEvents.mat';
load(path2load)

% Lowbeta

Indi1 = 1:2:10;
Indi2 = 2:2:10;

datasel1(1) = squeeze(nanmean(nanmean(nanmean(Peak.ERD.lowbeta.STD(Indi1,:,:)))));
datasel1(3) = squeeze(nanmean(nanmean(nanmean(Peak.ERS.lowbeta.STD(Indi1,:,:)))));

datasel2(1) = squeeze(nanmean(nanmean(nanmean(Peak.ERD.lowbeta.STD(Indi2,:,:)))));
datasel2(3) = squeeze(nanmean(nanmean(nanmean(Peak.ERS.lowbeta.STD(Indi2,:,:)))));

datasel1 = [0 0 datasel1 0 0];
datasel2 = [0 0 datasel2 0 0];

% Stem Plots

Colors = {'b', 'r'};

figure('color', 'white');
set(gcf, 'Position', [680 558 370 420]);

stem(datasel1(1:2:end), 'filled', Colors{1}, 'linewidth', 2); hold on
stem((1:round(length(datasel2)/2))+0.2,datasel2(1:2:end), 'filled', Colors{2}, 'linewidth', 2);

% Take out 0 points
% Beginning
stem(datasel1(1), 'filled', '-w', 'linewidth', 2);
stem(1+0.2,datasel2(1), 'filled', '-w', 'linewidth', 2);
% End
stem(4,datasel1(end), 'filled', '-w', 'linewidth', 2);
stem(4+0.2,datasel2(end), 'filled', '-w', 'linewidth', 2);

% Last touch
box off
a = gca; a.XColor = [1 1 1];
a.FontSize = 14;
a.FontWeight = 'bold';
a.YLabel.String = 'Norm %change';
hline = refline([0,0]);
hline.Color = 'k';
h = legend;
h.String{1} = 'STDs';
h.String{2} = 'STDw';
h.String(3:end) = [];
legend boxoff


% Violin plots

% Pre-stim

datasel1 = squeeze(nanmean(nanmean(Peak.ERD.lowbeta.STD(Indi1,:,:))));
datasel2 = squeeze(nanmean(nanmean(Peak.ERD.lowbeta.STD(Indi2,:,:))));

figure('color', 'white')
set(gcf, 'Position', [680 687 560 291]);

distributionPlot([datasel1,datasel2], 'color', {'b', 'r'})
title('Pre-')
set(gca, 'fontsize', 14, 'fontweight', 'bold')

h = gca; h = h.Children(3); % first event
h2 = gca; h2 = h2.Children(4); % second event

leg = legend([h, h2],{'STDw', 'STDs'});
leg.AutoUpdate = 'off';
leg.Box = 'off';


% Post-stim

datasel1 = squeeze(nanmean(nanmean(Peak.ERS.lowbeta.STD(Indi1,:,:))));
datasel2 = squeeze(nanmean(nanmean(Peak.ERS.lowbeta.STD(Indi2,:,:))));

figure('color', 'white')
set(gcf, 'Position', [680 687 560 291]);

distributionPlot([datasel1,datasel2], 'color', {'b', 'r'})
title('Post-')
set(gca, 'fontsize', 14, 'fontweight', 'bold')

h = gca; h = h.Children(3); % first event
h2 = gca; h2 = h2.Children(4); % second event

leg = legend([h, h2],{'STDw', 'STDs'});
leg.AutoUpdate = 'off';
leg.Box = 'off';


%% FFT plot

path2load = 'E:\DATA\Humans - Beatit\DATAproc\FFT\GA FFT.mat';
load(path2load)

FqOI = [1 4];
FqInd = find(ISO.xaxis >= FqOI(1) & ISO.xaxis <= FqOI(2));

figure('Color', 'White'); hold on, box off
set(gcf,  'Position', [680 558 419 420], 'name', 'GA FFT');

shadedErrorBar(ISO.xaxis(FqInd), ISO.nanmean2(FqInd), ISO.SE2(FqInd), '-b'); %'lineProps',
% title(['GA, FFT ISO fqcorr'])

a = gca;
a.XAxis.TickValues = 1 : .5 : FqOI(2);
a.FontSize = 14;
a.FontWeight = 'bold';

xlabel('\fontsize{14} Frequency (Hz)', 'fontweight', 'bold')
ylabel('\fontsize{14} Normalized Pow', 'fontweight', 'bold')

dir2save = 'C:\Users\p70068941\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Figures';
saveFigures_AC(dir2save)

%% DEV ERP for Binary Beat

load("E:\DATA\Humans - Beatit\DATAproc\ERP\meancorr\Results\DEV_Beat.mat")

figure; fullScreenFig

EVENpos = squeeze(nanmean(TwoBinary.ERP(:,[1,3], :),2));
EVENposSE = std(squeeze(nanmean(TwoBinary.ERP(:,[1,3], :),2)), 'omitnan')/sqrt(size(TwoBinary.ERP,1));

ODDpos = squeeze(nanmean(TwoBinary.ERP(:,[1,3]+1, :),2));
ODDposSE = std(squeeze(nanmean(TwoBinary.ERP(:,[1,3]+1, :),2)), 'omitnan')/sqrt(size(TwoBinary.ERP,1));

subplot(1,3,1:2)

shadedErrorBar(TwoBinary.Time, nanmean(EVENpos), EVENposSE, '-b'); hold on
shadedErrorBar(TwoBinary.Time, nanmean(ODDpos), ODDposSE, '-r'); hold on

AddLines_AC; AdjAxes_AC; box off

xlabel('\fontsize{14} Time (s)', 'fontweight', 'bold')
ylabel('\fontsize{14} Amplitude uV', 'fontweight', 'bold')

% Legend

h = gca; idx = length(h.Children);
h = h.Children(idx(1)); % first event
h2 = gca; h2 = h2.Children(idx-4); % second event

leg = legend([h, h2], {'Odd', 'Even'});
leg.AutoUpdate = 'off';
leg.Box = 'off';

title('Average Binary')

% Extract Amps

subplot(1,3,3); hold on

TimeOI = [.28 .31];
TimeInd = find(TwoBinary.Time >= TimeOI(1) & TwoBinary.Time <= TimeOI(2));

EVENposAmp = nanmean(EVENpos(:,TimeInd),2);
ODDposAmp = nanmean(ODDpos(:,TimeInd),2);

distributionPlot([EVENposAmp ODDposAmp], 'color', {'b', 'r'});
legend('Odd', 'Even'); legend boxoff; box off
leg.AutoUpdate = 'off';
title('P3 Amp')
AdjAxes_AC
YLIM = get(gca, 'ylim'); YLIM = YLIM(2);

% Stats

clear p 
[H, p, CI, Stats] = ttest(EVENposAmp, ODDposAmp);

% Plot stats
Stars = nan([1,2]);
if p < .05
    Stars(1) = YLIM*.8;
    plot(1:2, Stars, '*k');
    set(legend, 'String', {'Odd','Even','p <.05'});
end
set(gca, 'ylim', [-10 10])

% SingleSubj plots

n = 1;

for ss = 1:2:size(TwoBinary.ERP,1)
    
    EVENpos = squeeze(nanmean(TwoBinary.ERP(ss:ss+1,[1,3], :),2));
    EVENposSE = std(squeeze(nanmean(TwoBinary.ERP(:,[1,3], :),2)), 'omitnan')/sqrt(size(TwoBinary.ERP,1));
    ODDpos = squeeze(nanmean(TwoBinary.ERP(ss:ss+1,[1,3]+1, :),2));
    ODDposSE = std(squeeze(nanmean(TwoBinary.ERP(:,[1,3]+1, :),2)), 'omitnan')/sqrt(size(TwoBinary.ERP,1));
   
    % Time-course
    figure; fullScreenFig
    shadedErrorBar(TwoBinary.Time, nanmean(EVENpos), EVENposSE, '-b'); hold on
    shadedErrorBar(TwoBinary.Time, nanmean(ODDpos), ODDposSE, '-r'); hold on
    
    AddLines_AC; AdjAxes_AC; box off; title(sprintf('Subj %d', n));
    
    xlabel('\fontsize{14} Time (s)', 'fontweight', 'bold')
    ylabel('\fontsize{14} Amplitude uV', 'fontweight', 'bold')

    % Extract Amps
    
    figure(2); fullScreenFig; subplot(4,5,n); hold on; 
   
    TimeOI = [.28 .31];
    TimeInd = find(TwoBinary.Time >= TimeOI(1) & TwoBinary.Time <= TimeOI(2));
    
    EVENposAmp = nanmean(EVENpos(:,TimeInd),2);
    ODDposAmp = nanmean(ODDpos(:,TimeInd),2);
    
    distributionPlot([EVENposAmp ODDposAmp], 'color', {'b', 'r'});
%     legend('Odd', 'Even'); legend boxoff; box off
    title(sprintf('Subj %d', n)); n = n+1;
    AdjAxes_AC
    YLIM = get(gca, 'ylim'); YLIM = YLIM(2);
    
    % Stats
    
    clear p; alpha = .05;
    [H, p, CI, Stats] = ttest(EVENposAmp, ODDposAmp);
    
    % Plot stats
    Stars = nan([1,2]);
    if p < alpha
        Stars(2) = YLIM*.8;
        h = plot(1:2, Stars, '*k');
       legend(h,'p <.05');
    end
%     set(gca, 'ylim', [-10 10])
end
    
    

%%

loadName = "E:\DATA\Humans - Beatit\DATAproc\TFR\meancorr\Results\lowbeta RegrModels.mat";

% Plot modelling info

Labels = {'Binary', 'Triplet', 'CombBeat', 'Other'};

load(loadName)
NSeq = 96*2;
Info2Display = (Info2Display / NSeq) *100;

figure('color', 'white');
Colors = {'b', 'c', 'g', 'm'};
distributionPlot(Info2Display, 'color', Colors); box off
set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1])
ylabel('% trials')
legend(Labels, 'location', 'northwestoutside'); legend boxoff

% Save Excel
saveName = replace(loadName, 'lowbeta RegrModels.mat', 'BeatModels.xls');
writetable(table2array(round(single(Info2Display),1), 'VariableNames', Labels), saveName)


