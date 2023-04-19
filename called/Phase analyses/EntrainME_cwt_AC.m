

function EntrainME_cwt_AC(sess)

if nargin < 1
    sess= input('Session number?');
end

% Prepare data

datadir = getSessDir(sess);
outdir = fullfile(datadir.sess, 'DATAproc'); %, 'TFR');
outdir1 = fullfile(outdir, 'FFT');
outdir2 = fullfile(outdir, 'Entrainment');

if ~exist(outdir2)
    mkdir(outdir2)
end
dir2save = 'C:\Users\p70068941\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Figures\Phase synch';

files = dir(fullfile(outdir1, 'ep4FFT *.mat'));
if isempty(files)
    error('Call Data4FFT first')
end

% Set up parameters

Plot = 0; % 1 = yes, plot

ISO = 0.65; %ms
Stimfreq = round(1/ISO, 2); %Hz
FqOI = [Stimfreq-.5 Stimfreq+.5];
datatype = 1; % 1 = complex exponentials; 0 = real numbers (but here we need 1)

FqLab = 'StimFreq';

% FqOI = [12 25]; % for extraction
% FqLab = 'Beta';

load(fullfile(outdir1, files(1).name))

% Ntones = data.Ntones; % max length of sequence
% StimRate = data.stimRate_Hz; %Hz
ISO = data.ISO;
% Samp = data.fsample;
Time = data.time{1};

PlotBins = 10;
Pos = 8; % max 8 pos

% TimeSeq = (1:Ntones)*ISO - ISO; % From first onset (0) til 13th tone onset
% TimeOI_Start = ISO * Samp; % TimeInd for the beginning of the seq

Ntrials = 192;

Pad = round(1/Stimfreq/10,2); % padding around tone onset - for phase extraction

% Channels of interest
chLabels = 'E:\DATA\Humans - Beatit\DATAproc\chs_OI.mat';
load(chLabels); %chs_OI

% Prepare file names

saveName2 = sprintf('Phase Stats %s.mat', FqLab);
saveName2 = fullfile(outdir2, saveName2);

saveName3 = sprintf('GroupData %s.mat', FqLab);
saveName3 = fullfile(outdir2, saveName3);

% Looping variables

n = 1; % will be used for subjs
items = 7; % per fig, horizontally
rows = 4; % per fig

subp = 1:items*rows; % organize figures, see below
subp = reshape(subp, [items,rows])'; % 4 sessions per figure (vertically)

% Pre-allocate variables to save

ANGLES1 = nan([datadir.NSubj/2, Ntrials, length(chs_OI), Pos]);
ANGLES2 = ANGLES1;
GA_MVL = nan([datadir.NSubj/2, Ntrials, length(chs_OI)]);
GA_MVLpost = GA_MVL; GA_MVLrand = GA_MVL;

if ~exist(saveName3)
    for ff = 1:2:length(files)
        
        fprintf('\nWorking on Subj %d\n', n)
        row = 1; % counters of row-elements in plots
        
        % Savename
        saveName = sprintf('%s Subj%d.mat', FqLab, n);
        saveName = fullfile(outdir2, saveName);
        
        if ~exist(saveName)
            
            % Load 2 sess for each subj
            fullfilename = fullfile(outdir1, files(ff).name);
            load(fullfilename);
            Events1 = data.events; % this will get lost in ft_preproc
            
            data1 = get_TFR_AC(fullfilename, outdir2, datatype, FqOI);
            
            % And sess2
            fullfilename = fullfile(outdir1, files(ff+1).name);
            load(fullfilename);
            Events2 = data.events; % this will get lost in ft_preproc
            
            data2 = get_TFR_AC(fullfilename, outdir2, datatype, FqOI);
            
            % Concatenate 2 sessions
      
            data1.fourierspctrm = cat(1, data1.fourierspctrm, data2.fourierspctrm);
            Events = [Events1, Events2];
            Time = data1.time;

            clear data2
            data = data1; clear data1

            % Select FC cluster
            % 1. get index for chsOI
            
            Labels = data.label;
            indch = [];
            
            for cc = 1:length(chs_OI)
                inds = find(strcmp(Labels, chs_OI{cc}));
                indch = [indch inds];
            end
            
            % 2. Move to 3D Data matrix
            % average over fqsOI and select only chsOI

            Data = squeeze(nanmean(data.fourierspctrm(:,indch,:,:),3));
            
            %     figure; plot(real(squeeze(Data(:,end,:))'))
            %     figure; plot(imag(squeeze(Data(:,end,:))'))
            %     figure; plot(angle(squeeze(Data(:,end,:))'))
            
            
            clear data; data = Data;
            save(saveName, 'data', 'Events', 'Time')
            
        else
            load(saveName)
        end
        
        if ff == 1
            DATA = nan([datadir.NSubj/2, Ntrials, length(chs_OI), size(data,3)]);
        end
        
        % Prepare plots
        
        % Phase angles
        angles = angle(data); % in radians
        angles1 = nan([size(data,1), size(data,2), Pos]); % pre-allocate pre-DEV
        angles2 = nan(size(angles1)); % post-DEV
        
        % Random angles
        Null = wrap(deg2rad(rand(size(angles1))*360));
        
        % Loop over the sequences to retrieve DEV pos
        
        for tr = 1:length(Events)
            
            Ntones = length(Events{tr});
            DEVpos = find(strcmp(Events{tr}, 'DEV'));
            
            TonesOI = 3:DEVpos-1; % STD
            TonesOI2 = DEVpos+1:Ntones-1; % post-DEV
            
            Onset = TonesOI*ISO - ISO; % pre-DEV
            Onset2 = TonesOI2*ISO - ISO; % post-DEV
            
            % TimeInd pre-DEV
            for on = 1:length(Onset) % starting from the third, actually
                tind = find(Time >= Onset(on)-Pad & Time <= Onset(on));
                angles1(tr,:,on) = circ_mean(angles(tr,:,tind), [],3); % pre-DEV
            end
            
            % TimeInd post-DEV
            for on = 1:length(Onset2)
                tind = find(Time >= Onset2(on)-Pad & Time <= Onset2(on));
                angles2(tr,:,on) = circ_mean(angles(tr,:,tind), [],3); % post-DEV
            end
            
            % Since we are looping, let's add sth
            % Mean vector length for angles1 and random angles
            for chs = 1:size(angles,2)
                MVL(tr,chs) = circ_r(squeeze(angles1(tr,chs,squeeze(abs(angles1(tr,chs,:))) > 0)));
                MVLpost(tr,chs) = circ_r(squeeze(angles2(tr,chs,squeeze(abs(angles2(tr,chs,:))) > 0)));
                MVLrand(tr,chs) = circ_r(squeeze(Null(tr,chs,squeeze(abs(Null(tr,chs,:))) > 0)));
            end
        end
        
        % Envelope if necessary
        
        if FqOI(1) > 3
            data = envelope(real(data),30, 'peak');
        else
            data = real(data);
        end
        
        % Entrainment plots
        
        for chs = 1:size(data,2)
            
            tmpangles1 = []; tmpangles1 = squeeze(angles1(:,chs,:));
            tmpangles2 = []; tmpangles2 = squeeze(angles2(:,chs,:));
            
            if row == 1
                figure('color', 'white');
                fullScreenFig
            end
            
            SEdata = (std(data(:,chs,:), 'omitnan'))/sqrt(size(data,1));
            
            subplot(rows,items,subp(row,1:4));
            fig = shadedErrorBar(Time, squeeze(nanmean(data(:,chs,:))), SEdata, 'lineprops', '-b');  %'lineprops', 'b'
            hold on; title(sprintf('Subj %d %s ch',n, chs_OI{chs})); box off; AdjAxes_AC
            
            % Add lines
            
            TonesOI = 2:7; % STD
            TonesOI2 = 8:13; % post-DEV
            
            Onset = TonesOI*ISO - ISO; % pre-DEV
            Onset2 = TonesOI2*ISO - ISO; % post-DEV
            
            ylimit = fig.mainLine.Parent.YLim;
            
            Onsets = [Onset' Onset'];
            Onsets2 = [Onset2' Onset2'];
            line(Onsets, ylimit, 'color','blue','LineStyle','--'); % STD
            line(Onsets2(1,:), ylimit, 'color','green','LineStyle','--'); % DEV
            line(Onsets2(2:end,:), ylimit, 'color','red','LineStyle','--'); % post-DEV
            
            % Polar plot
            
            subplot(rows,items,subp(row,5));
            polarhistogram(tmpangles1(abs(tmpangles1) > 0),PlotBins, 'facecolor', 'b');
            a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold';
            title('Pre-DEV', 'fontsize', 12); % AdjAxes_AC
            
            subplot(rows,items,subp(row,6));
            polarhistogram(tmpangles2(abs(tmpangles2) > 0),PlotBins, 'facecolor', 'r');
            a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold';
            title('Post-DEV', 'fontsize', 12); % AdjAxes_AC
            
            % Plot mean vector lengths
            
            Nbins = PlotBins; % or different?
            
            subplot(rows,items,subp(row,end));
            histogram(MVL(:,chs), Nbins, 'facecolor', 'b'); hold on;
            h1 = histogram(MVLrand(:,chs), Nbins, 'facecolor', 'c');
            histogram(MVLpost(:,chs), Nbins, 'facecolor', 'r'); box off;
            a = gca; a.XLim = [0 1]; a.XAxis.TickValues = [0 .5 1]; a.YLim = [0 100];
            AdjAxes_AC
            
            if row == 1
                %             legend('pre-DEV', 'rand', 'post-DEV', 'Location', 'northeast', 'fontsize', 9); legend boxoff
                legend(h1, 'rand', 'Location', 'northeast', 'fontsize', 9); legend boxoff
            end
            title(sprintf('Mean Vector Length\nTrial-based'), 'fontsize', 12)
            
            row = row+1;
            
            % Stats with circular statistics toolbox and permutation tests
            
            permutations = 1000;
            try, [p, obsdiff, effectsize, ~] = circ_permutationTest_AC(tmpangles1,tmpangles2, permutations);
                
                % Store into a cell struct
                Stats.(chs_OI{chs}){1,1} = 'pvals'; Stats.(chs_OI{chs}){1,2} = 'Obs Diff'; Stats.(chs_OI{chs}){1,3} = 'Eff Size';
                Stats.(chs_OI{chs}){n+1,1} = p; Stats.(chs_OI{chs}){n+1,2} = obsdiff; Stats.(chs_OI{chs}){n+1,3} = effectsize;
            end
            
            % Stats on angles1 compared to random distribution
            
            try,  [p, obsdiff, effectsize] = circ_permutationTest_AC(tmpangles1, squeeze(Null(:,chs)), permutations);
                
                % Store into a cell struct
                StatsNull.(chs_OI{chs}){1,1} = 'pvals'; StatsNull.(chs_OI{chs}){1,2} = 'Obs Diff'; StatsNull.(chs_OI{chs}){1,3} = 'Eff Size';
                StatsNull.(chs_OI{chs}){n+1,1} = p; StatsNull.(chs_OI{chs}){n+1,2} = obsdiff; StatsNull.(chs_OI{chs}){n+1,3} = effectsize;
            end
            
            % Stats on mean vector lengths
            [p, obsdiff, effectsize] = permutationTest(MVL(:,chs),MVLrand(:,chs), permutations);
            
            a = gca;
            if p < .001
                p2show = 'p < .001';
            else
                p2show = sprintf('p = %s', num2str(round(p,3)));
            end
            text(a.XLim(2)*.5, a.YLim(2)*.7, p2show, 'fontsize', 12, 'fontweight', 'bold', 'color', 'c');
            
            StatsMVL.(chs_OI{chs}){1,1} = 'pvals'; StatsMVL.(chs_OI{chs}){1,2} = 'Obs Diff'; StatsMVL.(chs_OI{chs}){1,3} = 'Eff Size';
            StatsMVL.(chs_OI{chs}){n+1,1} = p; StatsMVL.(chs_OI{chs}){n+1,2} = obsdiff; StatsMVL.(chs_OI{chs}){n+1,3} = effectsize;
            
            % And pre- vs post-DEV
            [p, obsdiff, effectsize] = permutationTest(MVL(:,chs),MVLpost(:,chs), permutations);
            
            a = gca;
            if p < .001
                p2show = 'p < .001';
            else
                p2show = sprintf('p = %s', num2str(round(p,3)));
            end
            text(a.XLim(2)*.5, a.YLim(2)*.6, p2show, 'fontsize', 12, 'fontweight', 'bold', 'color', 'r');
            
            StatsMVLpost.(chs_OI{chs}){1,1} = 'pvals'; StatsMVLpost.(chs_OI{chs}){1,2} = 'Obs Diff'; StatsMVLpost.(chs_OI{chs}){1,3} = 'Eff Size';
            StatsMVLpost.(chs_OI{chs}){n+1,1} = p; StatsMVLpost.(chs_OI{chs}){n+1,2} = obsdiff; StatsMVLpost.(chs_OI{chs}){n+1,3} = effectsize;
            
            if row > 4, row = 1; end
            
            % Store Group data
            
            DATA(n,:,chs,:) = squeeze(data(:,chs,:));
            ANGLES1(n,:,chs,:) = tmpangles1;
            ANGLES2(n,:,chs,:) = tmpangles2;
            GA_MVL(n,:,chs) = MVL(:,chs);
            GA_MVLpost(n,:,chs) = MVLpost(:,chs);
            GA_MVLrand(n,:,chs) = MVLrand(:,chs);
            
        end
        
        % Calculate MVL for group data
        
        if n == datadir.NSubj/2
            MVLavg = squeeze(nanmean(GA_MVL,3));
            MVLavgPost = squeeze(nanmean(GA_MVLpost,3));
            MVLavgRand = squeeze(nanmean(GA_MVLrand,3));
        end
        
        % Save figures
        dirtmp = fullfile(dir2save, sprintf('Subj%d', n));
        saveFigures_AC(dirtmp); close all
        % Next subj
        n = n+1;
    end
    
    % Save Group data and Stats
    
    save(saveName3, 'DATA', 'ANGLES1', 'ANGLES2', 'GA_MVL', 'GA_MVLpost', 'GA_MVLrand', 'MVLavg', 'MVLavgPost', 'MVLavgRand', '-v7.3');
    %     save(saveName2, 'Stats', 'StatsNull', 'StatsMVL', 'StatsMVLpost');
    save(saveName2, 'StatsMVL', 'StatsMVLpost');
else
    load(saveName2)
    load(saveName3)
end


if Plot
    
    % Preapre GA plots
    for chs = 1:size(DATA,3)
        
        tmpangles1 = squeeze(ANGLES1(:,:,chs,:));
        tmpangles2 = squeeze(ANGLES2(:,:,chs,:));
        
        figure('color', 'white');
        set(gcf, 'Position', [11 708 1904 287], 'Name', 'GA');
        
        % Entrainment plots
        
        MEANdata = squeeze(nanmean(nanmean(squeeze(DATA(:,:,chs,:)))));
        SEdata = (std(nanmean(squeeze(DATA(:,:,chs,:))), 'omitnan'))/sqrt(size(DATA,1));
        
        subplot(1,items,1:4);
        fig = shadedErrorBar(Time, MEANdata, SEdata, '-b');  %'lineprops',
        hold on; box off
        AdjAxes_AC
        
        % Add lines
        
        ylimit = fig.mainLine.Parent.YLim;
        
        TonesOI = 2:7; % STD
        TonesOI2 = 8:13; % post-DEV
        
        Onset = TonesOI*ISO - ISO; % pre-DEV
        Onset2 = TonesOI2*ISO - ISO; % post-DEV
        
        Onsets = [Onset' Onset'];
        Onsets2 = [Onset2' Onset2'];
        line(Onsets, ylimit, 'color','blue','LineStyle','--'); % STD
        line(Onsets2(1,:), ylimit, 'color','green','LineStyle','--'); % DEV
        line(Onsets2(2:end,:), ylimit, 'color','red','LineStyle','--'); % post-DEV
        
        % Polar plot
        
        PlotBins = 7;
        
        subplot(1,items,5);
        polarhistogram(tmpangles1(abs(tmpangles1) > 0),PlotBins, 'FaceColor','b');
        a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold'; % AdjAxes_AC
        title(sprintf('Pre-DEV %s ch', chs_OI{chs}), 'fontsize', 12)
        
        subplot(1,items,6);
        polarhistogram(tmpangles2(abs(tmpangles2) > 0),PlotBins, 'FaceColor','r');
        a = gca; a.RLim = a.RLim * 1.2; a.FontWeight = 'bold'; % AdjAxes_AC
        title(sprintf('Post-DEV %s ch', chs_OI{chs}), 'fontsize', 12)
        
        % MVL plot
        
        nBins = 5;
        
        tmpMVLavg = squeeze(nanmean(GA_MVL(:,:,chs)));
        tmpMVLavgPost = squeeze(nanmean(GA_MVLpost(:,:,chs)));
        tmpMVLrand = MVLrand(:,chs);
        
        subplot(1,items,7); hold on
        histogram(tmpMVLrand, nBins, 'facecolor', 'c');
        h1 = histogram(tmpMVLavgPost, nBins, 'facecolor', 'r');
        histogram(tmpMVLavg, nBins, 'facecolor', 'b'); box off;
        a = gca; a.XAxis.TickValues = [0 .5 1]; a.YLim = a.YLim*1.3;
        AdjAxes_AC
        
        %     legend('rand', 'post-DEV', 'pre-DEV', 'location', 'northeast', 'fontsize', 9); legend boxoff
        legend(h1, 'rand', 'location', 'northeast', 'fontsize', 9); legend boxoff
        title(sprintf('Mean Vector Length\nTrial-based'), 'fontsize', 12)
        
        % Stats
        
        [p, obsdiff, effectsize, p2adjust] = permutationTest(tmpMVLavg,tmpMVLrand, permutations);
        
        if p < .001
            p2show = 'p < .001';
        else
            p2show = sprintf('p = %s', num2str(round(p,3)));
        end
        
        a = gca;
        text(a.XLim(2)*.5, a.YLim(2)*.7, p2show, 'fontsize', 12, 'fontweight', 'bold', 'color', 'c');
        
        % And pre-post DEV
        [p, obsdiff, effectsize, p2adjust] = permutationTest(tmpMVLavg,tmpMVLavgPost, permutations);
        
        if p < .001
            p2show = 'p < .001';
        else
            p2show = sprintf('p = %s', num2str(round(p,3)));
        end
        
        a = gca;
        text(a.XLim(2)*.5, a.YLim(2)*.65, p2show, 'fontsize', 12, 'fontweight', 'bold', 'color', 'r');
        
        % Save figures
        saveFigures_AC(fullfile(dir2save, sprintf('GA %s ch', chs_OI{chs}))); close all
        
        % And Stats tabs
        
        %     writecell(Stats, fullfile(dir2save, 'Stats Pre-Post-DEV.xlsx'));
        %     writecell(StatsNull, fullfile(dir2save, 'Stats Pre-Null.xlsx'));
        writecell(StatsMVL.(chs_OI{chs}), fullfile(dir2save, 'Stats MVL pre-DEV - Null.xlsx'), 'Sheet', chs);
        writecell(StatsMVLpost.(chs_OI{chs}), fullfile(dir2save, 'Stats MVL pre-post DEV.xlsx'), 'Sheet', chs);
        
    end
end

% MVLs over subjs and chs

figure; set(gcf, 'name', 'MVL Pooling across subj and chs', 'color', 'white', 'Position', [448 459 304 289]);

nBins = 7;

tmpMVLavg = reshape(MVLavg, [numel(MVLavg),1]); 
% tmpMVLavgPost = reshape(MVLavgPost, numel(MVLavgPost));
tmpMVLrand = reshape(MVLavgRand, [numel(MVLavgRand),1]);

histogram(tmpMVLavg, nBins, 'facecolor', 'b'); hold on; 
histogram(tmpMVLrand, nBins, 'facecolor', 'c'); box off;
% h1 = histogram(tmpMVLavgPost, nBins, 'facecolor', 'r');
a = gca; a.XAxis.TickValues = [0 .5 1]; a.YLim = a.YLim*1.3;
AdjAxes_AC

legend('pre-DEV', 'rand', 'location', 'northeast', 'fontsize', 9); legend boxoff
% legend(h1, 'MVLpost', 'location', 'northeast', 'fontsize', 9); legend boxoff
title(sprintf('Mean Vector Length\nTrial-based'), 'fontsize', 12)

% Stats

permutations = 1000;
[p, obsdiff, effectsize, p2adjust] = permutationTest(tmpMVLavg,tmpMVLrand, permutations);

if p < .001
    p2show = 'p < .001';
else
    p2show = sprintf('p = %s', num2str(round(p,3)));
end

a = gca;
text(a.XLim(2)*.6, a.YLim(2)*.7, p2show, 'fontsize', 12, 'fontweight', 'bold', 'color', 'k');

% Save
saveFigures_AC(dir2save); close all

% More figures

More_PhasePlots_AC(saveName3) % ssubj, trial-level
More_PhasePlots_ChsOI_AC(saveName3) % same as above, but keeping chs OI
PhasePlots_PoolingChsOI_AC(saveName3) % pooling chs OI



