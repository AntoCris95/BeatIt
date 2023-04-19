

% Created in November 2021
% Written by Antonio Criscuolo


function Plots_ERP_BinaryBeat(sess, BSLcorr)

% Set up

if nargin < 1
    sess = input('Session number?');
    BSLcorr= input('Use BSL-corrected data? 0= No, 1= Yes');
end

if BSLcorr == 0, disp('Using non-bsl corrected data'),
else disp('Using bsl-corrected data'), end


% Set up folders

datadir = getSessDir(sess);

if BSLcorr == 1
    outdir= fullfile(datadir.output, 'ERP', 'BSLcorr');
    lab = 'bslcorr';
else
    outdir= fullfile(datadir.output, 'ERP', 'meancorr');
    lab = 'meancorr';
end
load(fullfile(datadir.main, 'ERPheader.mat'));
NSubj = datadir.NSubj;
Fqs= 1;

datadir = fullfile(datadir.output, 'ERP');
saveName = fullfile(outdir, 'Results', 'new_DEV_Beat.mat');

if ~strcmp(getenv('COMPUTERNAME'), 'ACRISCUOLO')
    outdir2 = fullfile('C:\Users\p70068941\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\New Figures', 'DEVBeat');
else
    outdir2 = fullfile('C:\Users\anton\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\New Figures', 'DEVBeat');
end
if ~exist(outdir2)
    mkdir(outdir2)
end

% We will focus on DEV tones only

Event = 'DEV';
Labels= {'DEVs', 'DEVw'};
EventCodes = {'82', '92','102','112'};

% Define Time-windows of interest

Time = linspace(header.xstart, header.xend, header.datasize(6));

% ERP
TOF= [-.1 0.35];
TimeInd = find(Time >= TOF(1) & Time <= TOF(2));
TimeERP= Time(TimeInd);

% P50
TOF1= [0.03 0.07];
TimeInd1= find(Time >= TOF1(1) & Time <= TOF1(2));
Time1= Time(TimeInd1);

% N100
TOF2= [0.09 0.12];
TimeInd2= find(Time >= TOF2(1) & Time <= TOF2(2));
Time2= Time(TimeInd2);

% P200
TOF3= [.18 0.23];
TimeInd3= find(Time >= TOF3(1) & Time <= TOF3(2));
Time3= Time(TimeInd3);

% P300
TOF4= [.28 0.31];
TimeInd4= find(Time >= TOF4(1) & Time <= TOF4(2));
Time4= Time(TimeInd4);


%% Plot data accordint to Beat categorization

Beats = {'BinaryBeat', 'TripleBeat', 'CombBeat', 'OtherBeat'};
GroupFile = 'E:\DATA\Humans - Beatit\DATAproc\TFR\meancorr\Results\lowbeta RegrModels.mat';

files1 = dir(fullfile(datadir, 'ds *.mat'));

s = 1;

if exist(GroupFile)
    if ~exist(saveName)
        
        load(GroupFile); % Beat models
        
        for ss = 1:2:NSubj
            
            SubName = sprintf('Subj%d', s);
            
            % Determine SequenceN for DEV tones
            
            load(fullfile(datadir, files1(ss).name));
            DEV1 = ERP_DEVPos_AC(data.events);
            
            load(fullfile(datadir, files1(ss+1).name));
            DEV2 = ERP_DEVPos_AC(data.events);
            
            DEV.POS1 = [DEV1.POS1 DEV2.POS1*2];
            DEV.POS2 = [DEV1.POS2 DEV2.POS2*2];
            DEV.POS3 = [DEV1.POS3 DEV2.POS3*2];
            DEV.POS4 = [DEV1.POS4 DEV2.POS4*2];
            
            clear DEV1 DEV2
            
            for ee = 1:length(EventCodes) % 4 DEV positions
                
                % get field names from DEV structure
                FieldName = fieldnames(DEV);
                FieldName = FieldName{ee};
                
                % Load-in ERP data and concatenate 2 sessions
                files = dir(fullfile(outdir, sprintf('%s %s *.mat',lab,EventCodes{ee})));
                
                load(fullfile(outdir, files(ss).name));
                tempdata = data;
                
                load(fullfile(outdir, files(ss+1).name));
                tempdata = vertcat(tempdata,data);
                
                % Select trials with Binary beat
                TrlOI = ismember(DEV.(FieldName),BinaryBeat.(SubName).seqN);
                
                % Get ERP amps
                
                DATA{ee}.(SubName).ERP(:,1:length(TimeInd)) = squeeze(tempdata(TrlOI,end,:,:,:,TimeInd));
                
                % Mean amp around time of interest
                DATA{ee}.(SubName).P50 = squeeze(nanmean(tempdata(TrlOI,end,:,:,:,TimeInd1),6));
                DATA{ee}.(SubName).N100 = squeeze(nanmean(tempdata(TrlOI,end,:,:,:,TimeInd2),6));
                DATA{ee}.(SubName).P200 = squeeze(nanmean(tempdata(TrlOI,end,:,:,:,TimeInd3),6));
                DATA{ee}.(SubName).P300 = squeeze(nanmean(tempdata(TrlOI,end,:,:,:,TimeInd4),6));
                
                % Peak amp in time of interest
%                 DATA{ee}.(SubName).P50 = squeeze(max(tempdata(TrlOI,end,:,:,:,TimeInd1),[],6));
%                 DATA{ee}.(SubName).N100 = squeeze(min(tempdata(TrlOI,end,:,:,:,TimeInd2),[],6));
%                 DATA{ee}.(SubName).P200 = squeeze(max(tempdata(TrlOI,end,:,:,:,TimeInd3),[],6));
%                 DATA{ee}.(SubName).P300 = squeeze(max(tempdata(TrlOI,end,:,:,:,TimeInd4),[],6));
                
                DATA{ee}.(SubName).Time = TimeERP;
                
                % Select trials with Binary beat or Inverse binary
                
                TrlOI1 = find(BinaryBeat.(SubName).Bcoeff > 0); % binary
                TrlOI2 = find(BinaryBeat.(SubName).Bcoeff < 0); % Inv binary
                
                TrlOI1 = ismember(DEV.(FieldName),BinaryBeat.(SubName).seqN(TrlOI1));
                TrlOI2 = ismember(DEV.(FieldName),BinaryBeat.(SubName).seqN(TrlOI2));
                
                NewOrder = [2,1,4,3]; % to combine pure and inverse binary
                
                % Get ERP amps
                try
                    DATA{ee}.(SubName).ERP_PureBinary(:,1:length(TimeInd)) = squeeze(tempdata(TrlOI1,end,:,:,:,TimeInd));
                catch
                    DATA{ee}.(SubName).ERP_PureBinary(:,1:length(TimeInd)) = squeeze(tempdata(TrlOI1,end,:,:,:,TimeInd))';
                end
                try
                    DATA{ee}.(SubName).ERP_InvBinary(:,1:length(TimeInd)) = squeeze(tempdata(TrlOI2,end,:,:,:,TimeInd));
                catch
                    DATA{ee}.(SubName).ERP_InvBinary(:,1:length(TimeInd)) = squeeze(tempdata(TrlOI2,end,:,:,:,TimeInd))';
                end
                
                DATA{ee}.(SubName).P50_PureBinary(1:length(find(TrlOI1)),1) = squeeze(nanmean(tempdata(TrlOI1,end,:,:,:,TimeInd1),6));
                DATA{ee}.(SubName).N100_PureBinary(1:length(find(TrlOI1)),1) = squeeze(nanmean(tempdata(TrlOI1,end,:,:,:,TimeInd2),6));
                DATA{ee}.(SubName).P200_PureBinary(1:length(find(TrlOI1)),1) = squeeze(nanmean(tempdata(TrlOI1,end,:,:,:,TimeInd3),6));
                DATA{ee}.(SubName).P300_PureBinary(1:length(find(TrlOI1)),1) = squeeze(nanmean(tempdata(TrlOI1,end,:,:,:,TimeInd4),6));
                
                DATA{ee}.(SubName).P50_InvBinary(1:length(find(TrlOI2)),1) = squeeze(nanmean(tempdata(TrlOI2,end,:,:,:,TimeInd1),6));
                DATA{ee}.(SubName).N100_InvBinary(1:length(find(TrlOI2)),1) = squeeze(nanmean(tempdata(TrlOI2,end,:,:,:,TimeInd2),6));
                DATA{ee}.(SubName).P200_InvBinary(1:length(find(TrlOI2)),1) = squeeze(nanmean(tempdata(TrlOI2,end,:,:,:,TimeInd3),6));
                DATA{ee}.(SubName).P300_InvBinary(1:length(find(TrlOI2)),1) = squeeze(nanmean(tempdata(TrlOI2,end,:,:,:,TimeInd4),6));
            end
            
            
            % Now prepare TwoBinary by combining Pure and Inverse binary
            
            for ee = 1:length(EventCodes) % 4 DEV positions
                
                DATA{ee}.(SubName).ERP_TwoBinary = [DATA{ee}.(SubName).ERP_PureBinary; DATA{NewOrder(ee)}.(SubName).ERP_InvBinary];
                DATA{ee}.(SubName).P50_TwoBinary = [DATA{ee}.(SubName).P50_PureBinary; DATA{NewOrder(ee)}.(SubName).P50_InvBinary];
                DATA{ee}.(SubName).N100_TwoBinary = [DATA{ee}.(SubName).N100_PureBinary; DATA{NewOrder(ee)}.(SubName).N100_InvBinary];
                DATA{ee}.(SubName).P200_TwoBinary = [DATA{ee}.(SubName).P200_PureBinary; DATA{NewOrder(ee)}.(SubName).P200_InvBinary];
                DATA{ee}.(SubName).P300_TwoBinary = [DATA{ee}.(SubName).P300_PureBinary; DATA{NewOrder(ee)}.(SubName).P300_InvBinary];
                
            end
            
            s = s+1; % subjN
        end
        
        save(saveName, 'DATA')
        
        
        %% Preparing for plotting
        
        fprintf('\nPreapring SSubj plots')
        
        NSubj = length(fieldnames(DATA{1}));
        Colors = {'b', 'c', 'g', 'm'};
        Positions = {'8th', '9th', '10th', '11th'};
        ERPS = {'P50', 'N100', 'P200', 'P300'};
        
        figure('Color', 'White'); fullScreenFig; hold on % for full ERPs
        figure('Color', 'White'); fullScreenFig; hold on % for full ERPs, metric binary
        figure('Color', 'White'); fullScreenFig; hold on % for full ERPs, inverse binary
        figure('Color', 'White'); fullScreenFig; hold on % for full ERPs, two binary
        
        
        for ss = 1:NSubj
            
            SubName = sprintf('Subj%d', ss);
            Time = DATA{1}.(SubName).Time;
            
            % Plot ERPs
            for ee= 1:length(DATA) % 4 DEV positions: 8-11
                
                % ALL binary
                figure(1); subplot(4,5,ss); hold on
                set(gcf, 'Name', 'ALLData')
                plot(Time, nanmean(DATA{ee}.(SubName).ERP), 'linewidth', 2)
                title(SubName)
                set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xlim', [-.1 .35])
                legend(Positions); legend boxoff
                
                % Store GAs
                ALLData.ERP(ss,ee,:) = nanmean(DATA{ee}.(SubName).ERP);
                ALLData.P50(ss,ee) = nanmean(DATA{ee}.(SubName).P50);
                ALLData.N100(ss,ee) = nanmean(DATA{ee}.(SubName).N100);
                ALLData.P200(ss,ee) = nanmean(DATA{ee}.(SubName).P200);
                ALLData.P300(ss,ee) = nanmean(DATA{ee}.(SubName).P300);
                
                % Metric binary
                figure(2); subplot(4,5,ss); hold on
                set(gcf, 'Name', 'Pure Binary')
                plot(Time, nanmean(DATA{ee}.(SubName).ERP_PureBinary), 'linewidth', 2)
                title(SubName)
                set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xlim', [-.1 .35])
                legend(Positions); legend boxoff
                
                % Store GAs
                BinaryData.ERP(ss,ee,:) = nanmean(DATA{ee}.(SubName).ERP_PureBinary);
                BinaryData.P50(ss,ee) = nanmean(DATA{ee}.(SubName).P50_PureBinary);
                BinaryData.N100(ss,ee) = nanmean(DATA{ee}.(SubName).N100_PureBinary);
                BinaryData.P200(ss,ee) = nanmean(DATA{ee}.(SubName).P200_PureBinary);
                BinaryData.P300(ss,ee) = nanmean(DATA{ee}.(SubName).P300_PureBinary);
                
                % Inverse binary
                figure(3); subplot(4,5,ss); hold on
                set(gcf, 'Name', 'Inverse Binary')
                plot(Time, nanmean(DATA{ee}.(SubName).ERP_InvBinary), 'linewidth', 2)
                title(SubName)
                set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xlim', [-.1 .35])
                legend(Positions); legend boxoff
                
                % Store GAs
                InvBinaryData.ERP(ss,ee,:) = nanmean(DATA{ee}.(SubName).ERP_InvBinary);
                InvBinaryData.P50(ss,ee) = nanmean(DATA{ee}.(SubName).P50_InvBinary);
                InvBinaryData.N100(ss,ee) = nanmean(DATA{ee}.(SubName).N100_InvBinary);
                InvBinaryData.P200(ss,ee) = nanmean(DATA{ee}.(SubName).P200_InvBinary);
                InvBinaryData.P300(ss,ee) = nanmean(DATA{ee}.(SubName).P300_InvBinary);
                
                % Two binary
                figure(4); subplot(4,5,ss); hold on
                set(gcf, 'Name', 'Two Binary')
                plot(Time, nanmean(DATA{ee}.(SubName).ERP_PureBinary), 'linewidth', 2)
                title(SubName)
                set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xlim', [-.1 .35])
                legend(Positions); legend boxoff
                
                % Store GAs
                TwoBinaryData.ERP(ss,ee,:) = nanmean(DATA{ee}.(SubName).ERP_TwoBinary);
                TwoBinaryData.P50(ss,ee) = nanmean(DATA{ee}.(SubName).P50_TwoBinary);
                TwoBinaryData.N100(ss,ee) = nanmean(DATA{ee}.(SubName).N100_TwoBinary);
                TwoBinaryData.P200(ss,ee) = nanmean(DATA{ee}.(SubName).P200_TwoBinary);
                TwoBinaryData.P300(ss,ee) = nanmean(DATA{ee}.(SubName).P300_TwoBinary);
                
            end
            
            
            % Plot ODD/EVEN positions
            
            Labels = {'ALLData', 'BinaryData', 'InvBinaryData', 'TwoBinaryData'};
            
            figInd = [5:8];
            
            namefields = fieldnames(DATA{1}.Subj1);
            DataOI = namefields(startsWith(namefields, 'ERP'));
            
            for dd = 1:length(figInd)
                
                EVENpos = [];
                ODDpos = [];
                
                figure(figInd(dd)); fullScreenFig; subplot(4,5,ss); hold on
                set(gcf, 'color', [1 1 1], 'Name', Labels{dd})
                
                EVENpos(:,1:length(Time)) = DATA{1}.(SubName).(DataOI{dd});
                EVENpos = [EVENpos; DATA{3}.(SubName).(DataOI{dd})];
                
                ODDpos(:,1:length(Time)) = DATA{2}.(SubName).(DataOI{dd});
                ODDpos = [ODDpos; DATA{4}.(SubName).(DataOI{dd})];
                
                plot(Time, nanmean(EVENpos), 'linewidth', 2)
                plot(Time, nanmean(ODDpos), 'linewidth', 2)
                title(SubName)
                set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xlim', [-.1 .35])
                legend({'Even', 'Odd'}); legend boxoff
                
            end
        end
        
        % Save
        
        save(saveName,'DATA', 'ALLData', 'BinaryData', 'InvBinaryData', 'TwoBinaryData');
        
    else
        load(saveName)
    end
    
    % Plot extracted amplitudes
    
    namefields = fieldnames(DATA{1}.Subj1);
    NSubj = length(fieldnames(DATA{1}));
    Time = DATA{1}.Subj1.Time;
    
    Colors = {'b', 'c', 'g', 'm'};
    Positions = {'8th', '9th', '10th', '11th'};
    ERPS = {'P50', 'N100', 'P200', 'P300'};
    
    for ss = 1:NSubj
        
        SubName = sprintf('Subj%d', ss);
        
        figure('Color', 'White'); fullScreenFig; hold on
        set(gcf, 'name', SubName); N1 = get(gcf, 'number');
        
        figure('Color', 'White'); fullScreenFig; hold on
        set(gcf, 'name', SubName);  N2 = get(gcf, 'number');
        
        figpos = 1;
        figpos2 = 1;
        
        for dd = 1:4 % ALLdata up to two binary
            
            for pos = 1:length(ERPS) % P50 to P300
                
                % Plot 8-11th position
                
                txt = namefields(startsWith(namefields, ERPS{pos})); ERPLabel = txt{dd};
                
                figure(N1); subplot(4,4,figpos)
                violinplot_ERP_Binary(DATA, SubName, ERPLabel)
                title(ERPS{pos})
                set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1], 'ylim', [-30 30])
                
                figpos = figpos+1;
                
                % Plot extracted amplitudes for ODD/EVEN positions
                
                txt = namefields(startsWith(namefields, ERPS{pos})); ERPLabel = txt{dd};
                
                figure(N2); subplot(4,4,figpos2)
                violinplot_ERPsw_Binary(DATA, SubName, ERPLabel)
                title(ERPS{pos})
                set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1], 'ylim', [-30 30])
                
                figpos2 = figpos2+1;
            end
        end
    end
    
    
    
    %% Prepare GA plots
    
    fprintf('\nPreapring GA plots')
    DataOI = {'ALLData', 'BinaryData', 'InvBinaryData', 'TwoBinaryData'};
    
    for dd = 1:length(DataOI) % the above DataOI
        
        Data2Use = eval(DataOI{dd}); Data2Use = Data2Use.ERP;
        
        figure('Color', 'White'); fullScreenFig; hold on
        set(gcf, 'name', DataOI{dd})
        
        Colors = {'-b', '-c', '-g', '-m'};
        
        for ee = 1:size(Data2Use,2)
            
            Data = squeeze(Data2Use(:,ee,:));
            
            data1= [];
            data1.time= Time;
            data1.nanmean= nanmean(Data);
            data1.SE= (std(Data, 'omitnan'))/sqrt(NSubj); % maybe ntrials as well in the denominator?
            
            
            % Single plots
            
            shadedErrorBar(data1.time, data1.nanmean, data1.SE, 'lineprops', Colors{ee});
            hold on; box off;
            
            % Legend
            if ee == 4
                h = gca; h = h.Children(13); % first position
                h2 = gca; h2 = h2.Children(9); % second pos
                h3 = gca; h3 = h3.Children(5); % third pos
                h4 = gca; h4 = h4.Children(1); % fourth pos
                
                leg = legend([h, h2, h3, h4],Positions);
                leg.AutoUpdate = 'off';
                leg.Box = 'off';
                
                % Add lines
                plot1=gcf; plot1= plot1.CurrentAxes;
                plot1.FontSize = 14;
                plot1.FontWeight = 'bold';
                
                Minvals= [];
                Maxvals= [];
                Minvals= min(plot1.YLim);
                Maxvals= max(plot1.YLim);
                
                ylimit = [min(Minvals) max(Maxvals)];
                %[min(Minvals)+0.01*min(Minvals) max(Maxvals)+0.01*max(Maxvals)];
                ylim(ylimit);
                
                hline = refline([0,0]);
                hline.Color = 'k';
                line([0,0], ylimit, 'color','black','LineStyle','--');
                line([0.3,0.3], ylimit, 'color','black','LineStyle','--');
                
                xlabel('\fontsize{14} Time (s)', 'fontweight', 'bold')
                ylabel('\fontsize{14} Amplitude uV', 'fontweight', 'bold')
                
            end
        end
        
        
        GA plots for ODD/EVEN positions
        
        figure('Color', 'White'); fullScreenFig; hold on
        set(gcf, 'name', DataOI{dd});
        
        Colors = {'-b', '-r'};
        
        % EVEN
        Data1 = squeeze(nanmean(Data2Use(:,[1,3],:),2));
        
        data1= [];
        data1.time= Time;
        data1.nanmean= nanmean(Data1);
        data1.SE= (std(Data1, 'omitnan'))/sqrt(NSubj+2);
        
        % ODD
        Data2 = squeeze(nanmean(Data2Use(:,[2,4],:),2));
        
        data2= [];
        data2.time= Time;
        data2.nanmean= nanmean(Data2);
        data2.SE= (std(Data2, 'omitnan'))/sqrt(NSubj+2);
        
        
        % Single plots
        
        shadedErrorBar(data1.time, data1.nanmean, data1.SE, Colors{1});
        hold on; box off;
        shadedErrorBar(data1.time, data2.nanmean, data2.SE, Colors{2});
        
        % Legend
        h = gca; h = h.Children(5); % first event
        h2 = gca; h2 = h2.Children(1); % second event
        
        leg = legend([h, h2],{'Even', 'Odd'});
        leg.AutoUpdate = 'off';
        leg.Box = 'off';
        
        
        % Add lines
        plot1=gcf; plot1= plot1.CurrentAxes;
        plot1.FontSize = 14;
        plot1.FontWeight = 'bold';
        
        Minvals= [];
        Maxvals= [];
        Minvals= min(plot1.YLim);
        Maxvals= max(plot1.YLim);
        
        ylimit = [min(Minvals) max(Maxvals)];
        %[min(Minvals)+0.01*min(Minvals) max(Maxvals)+0.01*max(Maxvals)];
        ylim(ylimit);
        
        hline = refline([0,0]);
        hline.Color = 'k';
        line([0,0], ylimit, 'color','black','LineStyle','--');
        line([0.3,0.3], ylimit, 'color','black','LineStyle','--');
        
        xlabel('\fontsize{14} Time (s)', 'fontweight', 'bold')
        ylabel('\fontsize{14} Amplitude uV', 'fontweight', 'bold')
        
        % Statistical testing
        
        alpha = .05; %alpha/nperm
        permutations = 1000;
        
        [H,P,CI,STATS] = ttest(Data1, Data2);
        
        % Plot stats
        Stars = nan(size(P));
        Stars(find(P <= alpha)) = ylimit(2)*0.9;
        h3 = plot(data1.time, Stars, '*k');
        %         leg.String{3} = 'pAdj <.05';
        legend([h, h2, h3], {'Even', 'Odd', 'pAdj <.05'});
        
        
        % Plot extracted amplitudes
        
        figure('Color', 'White'); fullScreenFig; hold on % for P50, N100, P200 extracted amps
        set(gcf, 'name', DataOI{dd})
        Colors = {'b', 'c', 'g', 'm'};
        
        clear Data Data1 Data2
        
        Data2Use = eval(DataOI{dd});
        Data{1} = Data2Use.P50; Data{2} = Data2Use.N100; Data{3} = Data2Use.P200; Data{4} = Data2Use.P300;
        
        for pos = 1:4
            
            subplot(1,4,pos)
            distributionPlot(Data{pos}, 'color', Colors);
            legend(Positions); legend boxoff
            title(ERPS{pos})
            set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1], 'ylim', [-12 12])
            
        end
        
        
        % Plot extracted amplitudes for ODD/EVEN positions
        
        figure('Color', 'White'); fullScreenFig; hold on % for P50, N100, P200, P300 extracted amps
        set(gcf, 'name', DataOI{dd})
        Colors = {'b', 'r'};
        
        clear Data
        Data{1}(:,1) = nanmean(Data2Use.P50(:,[1,3]),2);
        Data{1}(:,2) = nanmean(Data2Use.P50(:,[1,3]+1),2);
        
        Data{2}(:,1) = nanmean(Data2Use.N100(:,[1,3]),2);
        Data{2}(:,2) = nanmean(Data2Use.N100(:,[1,3]+1),2);
        
        Data{3}(:,1) = nanmean(Data2Use.P200(:,[1,3]),2);
        Data{3}(:,2) = nanmean(Data2Use.P200(:,[1,3]+1),2);
        
        Data{4}(:,1) = nanmean(Data2Use.P300(:,[1,3]),2);
        Data{4}(:,2) = nanmean(Data2Use.P300(:,[1,3]+1),2);
        
        for pos = 1:4
            
            subplot(1,4,pos)
            distributionPlot(Data{pos}, 'color', Colors); hold on
            legend({'Even', 'Odd'}); legend boxoff
            title(ERPS{pos})
            set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1], 'ylim', [-12 12])
            
            % Stats
            
            [p, observeddifference, effectsize, p2Adjust] = permutationTest(Data{pos}(:,1), Data{pos}(:,2), permutations);
            
            % Plot stats
            Stars = nan([1,3]);
            if find(p < alpha)
                Stars(2) = 11;
                plot(1:3, Stars, '*k');
                %                  leg.String{3} = 'pAdj <.05';
                legend({'Even', 'Odd','pAdj <.05'});
            end
            
            clear p observeddifference effectsize p2Adjust Data1 data
            
        end
    end
    
    
    % Average Binary and Inv Binary
    
    fields2use = fieldnames(BinaryData);
    newOrd = [2,1,4,3]; % re-order even-odd positions
    
    TwoBinary.Time = TimeERP;
%     TwoBinary.ERP = [BinaryData.ERP; InvBinaryData.ERP(:,newOrd,:)];
%     
%     for ffs = 2:length(fields2use)
%         TwoBinary.(fields2use{ffs}) = [BinaryData.(fields2use{ffs}); InvBinaryData.(fields2use{ffs})(:,newOrd)];
%     end
    
    
    % Plot 8-11th positions
    
    figure; fullScreenFig; hold on; box off;
    Colors = {'b', 'c', 'g', 'm'};
    
    NSubj = size(TwoBinary.ERP,1);
    
    data2= [];
    data2.time= TimeERP;
    data2.nanmean= squeeze(nanmean(TwoBinary.ERP));
    data2.SE= squeeze((std(TwoBinary.ERP, 'omitnan')))/sqrt(NSubj);
    
    
    % Shaded plots
    
    for tr = 1:size(data2.nanmean,1)
        shadedErrorBar(data2.time, data2.nanmean(tr,:), data2.SE(tr,:), Colors{tr});
    end
    
    
    % Legend
    idx = tr*4; idx = 4:4:idx;
    h = gca; h = h.Children(idx(1)); % first event
    h2 = gca; h2 = h2.Children(idx(2)); % second event
    h3 = gca; h3 = h3.Children(idx(3)); % third event
    h4 = gca; h4 = h4.Children(idx(4)); % fourth event
    
    leg = legend([h, h2, h3, h4], {'8th', '9th', '10th', '11th'}); %{'Even', 'Odd'});
    leg.AutoUpdate = 'off';
    leg.Box = 'off';
    
    title('Average Binary')
    
    AddLines_AC
    AdjAxes_AC
    
    xlabel('\fontsize{14} Time (s)', 'fontweight', 'bold')
    ylabel('\fontsize{14} Amplitude uV', 'fontweight', 'bold')
    
    
    % Plot ODD/EVEN positions
    
    Colors = {'b', 'r'};
    
    ODD.nanmean = nanmean(data2.nanmean(1:2:end,:));
    ODD.SE = (std(squeeze(nanmean(TwoBinary.ERP(:,[1,3],:),2)), 'omitnan'))/sqrt(NSubj+2);
    
    EVEN.nanmean = nanmean(data2.nanmean(2:2:end,:));
    EVEN.SE = (std(squeeze(nanmean(TwoBinary.ERP(:,[1,3]+1,:),2)), 'omitnan'))/sqrt(NSubj+2);
    
    figure; fullScreenFig;
    shadedErrorBar(data2.time, ODD.nanmean, ODD.SE, Colors{1});
    shadedErrorBar(data2.time, EVEN.nanmean, EVEN.SE, Colors{2});
    
    
    % Legend
    idx = 2*4; idx = 2:4:idx;
    h = gca; h = h.Children(idx(1)); % first event
    h2 = gca; h2 = h2.Children(idx(2)); % second event
    
    leg = legend([h, h2], {'Odd', 'Even'});
    leg.AutoUpdate = 'off';
    leg.Box = 'off';
    
    title('Average Binary')
    
    AddLines_AC
    AdjAxes_AC
    %     set(gca, 'ylim', [-6 6]);
    
    xlabel('\fontsize{14} Time (s)', 'fontweight', 'bold')
    ylabel('\fontsize{14} Amplitude uV', 'fontweight', 'bold')
    
    
    % Run stats
    
    TimeTest = [.08 .2];
    TimeIndTest = find(data2.time >= TimeTest(1) & data2.time <= TimeTest(2));
    TimeTest = data2.time(TimeIndTest);
    
    [H,p,CI,Stats] = ttest(squeeze(nanmean(TwoBinary.ERP(:,[1,3],TimeIndTest),2)), squeeze(nanmean(TwoBinary.ERP(:,[1,3]+1,TimeIndTest),2)));
    
    % Plot stats
    Stars = nan([1,length(p)]);
    alpha = .05;
    if find(p < alpha)
        val = get(gca, 'ylim'); val = 0.8*val(2);
        Stars(find(p < alpha)) = val;
        h3 = plot(TimeTest, Stars, '*k');
        legend([h,h2,h3],{'Even', 'Odd','pAdj <.05'});
    end
    set(gcf, 'name', 'S-w real comparison')
    
    % Plot extracted amplitudes from 8-11th positions
    
    figure; fullScreenFig; hold on % for P50, N100, P200 extracted amps
    set(gcf, 'name', 'Average Binary')
    Colors = {'b', 'c', 'g', 'm'}; % Colors = {'b', 'r'};
    
    clear Data
    Data{1} = TwoBinary.P50; Data{2} = TwoBinary.N100; Data{3} = TwoBinary.P200; Data{4} = TwoBinary.P300;
    
    for pos = 1:4
        
        subplot(1,4,pos)
        distributionPlot(Data{pos}, 'color', Colors);
        legend(Positions); legend boxoff
        title(ERPS{pos})
        AdjAxes_AC
        set(gca, 'ylim', [-20 20])
        
    end
    
    % Plot extracted amplitudes for Even - Odd positions
    
    figure; fullScreenFig; hold on % for P50, N100, P200, P300 extracted amps
    set(gcf, 'name', 'Average Binary')
    Colors = {'b', 'r'};
    
    clear Data p
    
    for pos = 1:4
        
        Data{pos}(:,1) = vertcat(TwoBinary.(ERPS{pos})(:,1), TwoBinary.(ERPS{pos})(:,3));
        Data{pos}(:,2) = vertcat(TwoBinary.(ERPS{pos})(:,2), TwoBinary.(ERPS{pos})(:,4));
        
        subplot(1,4,pos)
        distributionPlot(Data{pos}, 'color', Colors);
        legend('Odd', 'Even'); legend boxoff
        title(ERPS{pos})
        AdjAxes_AC
        set(gca, 'ylim', [-20 20])
        
        % Stats
        
        p{pos} = nan([1,NSubj/2]);
        [H, p{pos}, CI, Stats] = ttest(Data{pos}(:,1), Data{pos}(:,2));
        
        % Plot stats
        Stars = nan([1,3]);
        if nanmean(p{pos}) < alpha
            Stars(2) = 15;
            plot(1:3, Stars, '*k');
            %                  leg.String{3} = 'pAdj <.05';
            legend({'Even', 'Odd','pAdj <.05'});
        end
        
        fprintf('\n%d Subj had a %s effect\n', sum(p{pos} < alpha), ERPS{pos});
    end
    
    save(saveName,'DATA', 'ALLData', 'BinaryData', 'InvBinaryData', 'TwoBinary');
    
    % Save all figures
    saveFigures_AC(outdir2); close all
    
    % Distribution of Binary and inverse binary
    
    outdir3 = 'C:\Users\p70068941\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Figures';
    
    % Select trials with Binary beat or Inverse binary
    
    SubNames = fieldnames(BinaryBeat);
    TotTrials = 192;
    TrlOI1 = [];
    TrlOI2 = [];
    
    for ss = 1:length(SubNames)
        
        TrlOI1(ss) = (length(find(BinaryBeat.(SubNames{ss}).Bcoeff > 0))/TotTrials)*100; % binary
        TrlOI2(ss) = (length(find(BinaryBeat.(SubNames{ss}).Bcoeff < 0))/TotTrials)*100; % Inv binary
        
    end
    
    figure; set(gcf, 'color', 'white', 'position', [403 246 325 420], 'name', 'BinaryBeats')
    distributionPlot([TrlOI1', TrlOI2'], 'color', {'b', 'c'})
    set(gca, 'fontsize', 14, 'fontweight', 'bold', 'ylim', [1 30])
    % legend('Binary', 'InvBinary', 'location', 'northwestoutside'); legend boxoff
    legend('S-w', 'w-S'); legend boxoff
    
    saveFigures_AC(outdir3); close
    
end
