
function PhasePlots_PoolingChsOI_AC(loadName)

load(loadName) 
% load(saveName3) % if debugging
PlotBins = 10;

% Channels of interest
chLabels = 'E:\DATA\Humans - Beatit\DATAproc\chs_OI.mat';
load(chLabels); %chs_OI

if ~strcmp(getenv('COMPUTERNAME'), 'ACRISCUOLO')
    outdir = fullfile('C:\Users\p70068941\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Figures\More Phase Synch\Pooling chsOI');
else
    outdir = fullfile('C:\Users\anton\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Figures\More Phase Synch\Pooling chsOI');
end

% Single-subj phase plots pooling over chs of interest
% Just a few samples from 192 trials (there would  be too many figures)

for ss = 1:size(ANGLES1,1)
    
    fprintf('\nWorking on Subj %d\n', ss)

    SubjName = sprintf('Subj%d', ss);
    outdir2 = fullfile(outdir, 'Raw Phase');
    
    tmptrls = sort(unique(randi(size(ANGLES1,2),[1,24]))); % 4 looping - will be used twice
    n = 1; % for subplotting
    for tr = tmptrls
        
        if n == 1
            figure; fullScreenFig;
            set(gcf, 'Name', SubjName)
        end
        
        angles1 = reshape(squeeze(ANGLES1(ss,tr,:,:)), [numel(squeeze(ANGLES1(ss,tr,:,:))), 1]);
        
        subplot(4,6,n);
        polarhistogram(angles1(abs(angles1) > 0),PlotBins, 'facecolor', 'b');
        a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold'; a.FontSize = 12;
        title(sprintf('Seq N%d', tr), 'fontsize', 12); 
        
        n = n+1;
        if n > 24
            n = 1;
            saveFigures_AC(outdir2); close all
        end
    end
    
    % Pool all trials and chs
    
    figure; set(gcf, 'name', ['Pooling ' SubjName], 'color', 'white', 'Position', [422 394 322 261]);

    angles1_tmp = reshape(ANGLES1(ss,:,:,:), [numel(ALLangles(ss,:,:,:)), 1]);
    
    polarhistogram(angles1_tmp(abs(angles1_tmp) > 0),PlotBins, 'facecolor', 'b');
    a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold'; a.FontSize = 12;
    title(['Pooling phase for ' SubjName], 'fontsize', 12);
    
    saveFigures_AC(outdir2); close all
    
    
    % Plot relative phase by pooling across chs and trials

    outdir2 = fullfile(outdir, 'Relative Phase');
    figure; set(gcf, 'name', ['Pooling ' SubjName], 'color', 'white', 'Position', [422 394 322 261]);

    angles1_tmp = reshape(squeeze(ALLangles(ss,:,:,:)), [numel(squeeze(ALLangles(ss,:,:,:))),1]);

    polarhistogram(angles1_tmp(abs(angles1_tmp) > 0),PlotBins, 'facecolor', 'b');
    a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold'; a.FontSize = 12;
    title(sprintf('Pooling phase for %s', SubjName), 'fontsize', 12);

    saveFigures_AC(outdir2); close all
end

% Finally, pool across subjs as well
% First for raw phase

outdir2 = fullfile(outdir, 'Raw Phase');
figure; set(gcf, 'name', 'Pooling all subjs', 'color', 'white', 'Position', [422 394 322 261]);

angles1_tmp = reshape(ANGLES1, [numel(ALLangles), 1]);

polarhistogram(angles1_tmp(abs(angles1_tmp) > 0),PlotBins, 'facecolor', 'b');
a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold'; a.FontSize = 12;
title('Pooling across subjs', 'fontsize', 12);

saveFigures_AC(outdir2); close all

% And for the relative phase

outdir2 = fullfile(outdir, 'Relative Phase');
figure; set(gcf, 'name', 'Pooling all subjs', 'color', 'white', 'Position', [422 394 322 261]);

angles1_tmp = reshape(ALLangles, [numel(ALLangles),1]);

polarhistogram(angles1_tmp(abs(angles1_tmp) > 0),PlotBins, 'facecolor', 'b');
a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold'; a.FontSize = 12;
title('Pooling across subjs', 'fontsize', 12);

saveFigures_AC(outdir2); close all
