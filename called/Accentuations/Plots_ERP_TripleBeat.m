

% Created in November 2021
% Written by Antonio Criscuolo


function Plots_ERP_TripleBeat(sess, BSLcorr)


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
saveName = fullfile(outdir, 'Results', 'DEV_TripleBeat.mat');

if ~strcmp(getenv('COMPUTERNAME'), 'ACRISCUOLO')
    outdir2 = fullfile('C:\Users\p70068941\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Figures', 'DEV TripleBeat');
else
    outdir2 = fullfile('C:\Users\anton\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Figures', 'DEV TripleBeat');
end
if ~exist(outdir2)
    mkdir(outdir2)
end

% We will focus on DEV tones only

Event = 'DEV';
Labels= {'DEVs', 'DEVw'};
EventCodes = {'82', '92','102','112'};
DEVpos = 8:11;
FirstDEV = 8; % first possible DEV pos
LastDEV = 11; % 11 tones up to the last DEV position

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

Beats = {'S-w-w', 'w-S-w', 'w-w-S'};
LNames = {'Triplet1', 'Triplet2', 'Triplet3'};

% Prepare accentuation patterns
Accents{1} = 1:3:LastDEV; Accents{1} = Accents{1}(Accents{1} >= FirstDEV);
Accents{2} = 2:3:LastDEV; Accents{2} = Accents{2}(Accents{2} >= FirstDEV);
Accents{3} = 3:3:LastDEV; Accents{3} = Accents{3}(Accents{3} >= FirstDEV);

% Files 2 load

GroupFile = 'E:\DATA\Humans - Beatit\DATAproc\TFR\meancorr\Results\TripleBeat_Model.mat';
files1 = dir(fullfile(datadir, 'ds *.mat'));

% Loop

DATA = []; % will store all data
s = 1; % looping var, keeps track of SubjN, despite the 2 sessions

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
        
        % Select trials with Triplet beat
        
        Triplets{1} = Triplet1.(SubName).seqN;
        Triplets{2} = Triplet2.(SubName).seqN;
        Triplets{3} = Triplet3.(SubName).seqN;
        
        DEVs = fieldnames(DEV);
        
        % Looping over 4 DEV positions
        
        AllData = []; % store data
        
        for ee = 1:length(EventCodes)
            
            % get field names from DEV structure
            FieldName = DEVs{ee};
            
            % Load-in ERP data and concatenate 2 sessions
            files = dir(fullfile(outdir, sprintf('%s %s *.mat',lab,EventCodes{ee})));
            
            load(fullfile(outdir, files(ss).name));
            AllData = [AllData; data];
            
            load(fullfile(outdir, files(ss+1).name));
            AllData = [AllData; data];
            
        end
        
        % Allocate DEV positions according to Accentuation types
        % Determine whether this is a S or w position
        % This will change according to Triplet type
        
        for ee = 1:length(EventCodes)
            
            S_Tr = []; % Strong trials
            w_Tr = []; % weak trials
            
            for accs = 1:length(Accents) % loop over 3 accents
                if any(Accents{accs} == DEVpos(ee)) % Accentuated, S
                    S_Tr = [S_Tr, Triplets{accs}];
                else % non-accentuated, w
                    w_Tr = [w_Tr, Triplets{accs}];
                end
            end
            
            % check N
            S_Tr = S_Tr(S_Tr <= size(AllData,1));
            w_Tr = w_Tr(w_Tr <= size(AllData,1));
            
            % Get ERPs
            
            DATA{ee}.(SubName).Time = TimeERP;
            DATA{ee}.(SubName).S_Tr = S_Tr;
            DATA{ee}.(SubName).w_Tr = w_Tr;
            
            DATA{ee}.(SubName).ERP_S = squeeze(AllData(S_Tr,end,:,:,:,TimeInd));
            DATA{ee}.(SubName).ERP_w = squeeze(AllData(w_Tr,end,:,:,:,TimeInd));
            
            % Peak amp in time of interest
            DATA{ee}.(SubName).P50_S = squeeze(max(AllData(S_Tr,end,:,:,:,TimeInd1),[],6));
            DATA{ee}.(SubName).N100_S = squeeze(min(AllData(S_Tr,end,:,:,:,TimeInd2),[],6));
            DATA{ee}.(SubName).P200_S = squeeze(max(AllData(S_Tr,end,:,:,:,TimeInd3),[],6));
            DATA{ee}.(SubName).P300_S = squeeze(max(AllData(S_Tr,end,:,:,:,TimeInd4),[],6));
            
            DATA{ee}.(SubName).P50_w = squeeze(max(AllData(w_Tr,end,:,:,:,TimeInd1),[],6));
            DATA{ee}.(SubName).N100_w = squeeze(min(AllData(w_Tr,end,:,:,:,TimeInd2),[],6));
            DATA{ee}.(SubName).P200_w = squeeze(max(AllData(w_Tr,end,:,:,:,TimeInd3),[],6));
            DATA{ee}.(SubName).P300_w = squeeze(max(AllData(w_Tr,end,:,:,:,TimeInd4),[],6));
            
        end
        s = s+1; % subjN
    end
    
    save(saveName, 'DATA')
    clear AllData
    
    %% Preparing GAs for plotting
    
    fprintf('\nPreapring GAs for plots')
    
    NSubj = length(fieldnames(DATA{1}));
    ERPS = {'P50', 'N100', 'P200', 'P300'};
    Beat_Pos = {'S', 'w'};
    
    ALLData = [];

    for ss = 1:NSubj

        SubName = sprintf('Subj%d', ss);
        ALLData.Time = DATA{1}.(SubName).Time;
        
        for ee = 1:length(DEVpos)
            for bps = 1:2 %ERP S,w
                
                FName = sprintf('ERP_%s', Beat_Pos{bps});
                ALLData.(FName)(ss,ee,:) = nanmean(DATA{ee}.(SubName).(FName));
                
                for erps = 1:4 % ERP comp, S,w, but also serves for 4 DEV pos
                    FName2 =  sprintf('%s_%s', ERPS{erps},Beat_Pos{bps});
                    ALLData.(FName2)(ss,ee) = nanmean(DATA{ee}.(SubName).(FName2));
                end
            end
        end
    end
    
    % Save
    
    save(saveName,'DATA', 'ALLData');

    else 
        load(saveName)
    end
    
    %% Plot S-w positions
    
    Data2plot.Time = ALLData.Time;
    Data2plot.ERP_S = squeeze(nanmean(nanmean(ALLData.ERP_S)));
    Data2plot.ERP_S_SE = squeeze((std(nanmean(ALLData.ERP_S,2))/sqrt(size(ALLData.ERP_S,1))));
    
    Data2plot.ERP_w = squeeze(nanmean(nanmean(ALLData.ERP_w)));
    Data2plot.ERP_w_SE = squeeze((std(nanmean(ALLData.ERP_w,2))/sqrt(size(ALLData.ERP_w,1))));
    
    figure; fullScreenFig
    set(gcf, 'name', 'GA S-w')
    
    shadedErrorBar(Data2plot.Time, Data2plot.ERP_S, Data2plot.ERP_S_SE, 'lineprops', 'b'); % 'lineprops'
    shadedErrorBar(Data2plot.Time, Data2plot.ERP_w, Data2plot.ERP_w_SE, 'lineprops', 'r');
    
    box off; title('DEV S-w')
    set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xlim', [-.1 .35])
    
    h = gca; h = h.Children(end); % first event
    h2 = gca; h2 = h2.Children(end-1); % second event, end-4
    a = gca;
    
    leg = legend([h, h2],{'S', 'w'});
    leg.AutoUpdate = 'off';
    leg.Box = 'off';
    
    [H,P,CI,STATS] = ttest(squeeze(nanmean(ALLData.ERP_S,2)), squeeze(nanmean(ALLData.ERP_w,2)));
    
    % Plot stats
    alpha = .05;
    Stars = nan(size(P));
    Stars(find(P <= alpha)) = a.YLim(2)*0.9;
    h3 = plot(ALLData.Time, Stars, '*k');
    %         leg.String{3} = 'pAdj <.05';
    legend([h, h2, h3], {'S', 'w', 'p <.05'});
    
    %% Plot extracted amplitudes
    
    namefields = fieldnames(DATA{1}.Subj1);
    NSubj = length(fieldnames(DATA{1}));
    
    Colors = {'b', 'c'};
    Positions = {'8th', '9th', '10th', '11th'};
    ERPS = {'P50', 'N100', 'P200', 'P300'};
    
    figure; fullScreenFig; set(gcf, 'name', 'ERP amps')
    
    for pos = 1:length(ERPS) % P50 to P300
        
        % Average across 8-11th position
        
        txt = sprintf('%s_%s', ERPS{pos}, 'S');
        txt2 = sprintf('%s_%s', ERPS{pos}, 'w');
        
        Data1 = nanmean(ALLData.(txt),2);
        Data2 = nanmean(ALLData.(txt2),2);
        
        subplot(1,4,pos); hold on
        distributionPlot([Data1, Data2], 'color', Colors);
        title(ERPS{pos});
        set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1])
        
        a= gca;
        
        % Statistical testing
        
        alpha = .05; %alpha/nperm
        permutations = 1000;
        [p, observeddifference, effectsize, p2Adjust] = permutationTest(Data1, Data2, permutations);
        
        [H,P,CI,STATS] = ttest(Data1, Data2);
        
        % Plot stats
        Stars = nan(size(P));
        if ~isnan(Stars), Stars(find(P <= alpha)) = a.YLim(2)*0.9;
            h3 = plot(Stars, '*k');
            %         leg.String{3} = 'pAdj <.05';
            legend({'S', 'w', 'p <.05'});
        else,  legend({'S', 'w'});  legend boxoff
        end
    end
   
    % Plot according to sequence position: 8-11th

    for pos = 1:length(ERPS) % P50 to P300

        figure; fullScreenFig; set(gcf, 'name', sprintf('%s amps', ERPS{pos}))

        for posi = 1:length(Positions) % 8-11th

            % Plot 8-11th position

            txt = sprintf('%s_%s', ERPS{pos}, 'S');
            txt2 = sprintf('%s_%s', ERPS{pos}, 'w');

            Data1 = ALLData.(txt)(:,posi);
            Data2 = ALLData.(txt2)(:,posi);

            subplot(1,4,posi); hold on
            distributionPlot([Data1, Data2], 'color', Colors);
            title(sprintf('%s %s pos', ERPS{pos}, Positions{posi}));
            set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1])

            a= gca;

            % Statistical testing

            alpha = .05; %alpha/nperm
            permutations = 1000;
            [p, observeddifference, effectsize, p2Adjust] = permutationTest(Data1, Data2, permutations);

            [H,P,CI,STATS] = ttest(Data1, Data2);

            % Plot stats
            Stars = nan(size(P));
            if ~isnan(Stars), Stars(find(P <= alpha)) = a.YLim(2)*0.9;
                h3 = plot(Stars, '*k');
                %         leg.String{3} = 'pAdj <.05';
                legend({'S', 'w', 'p <.05'});
            else,  legend({'S', 'w'});  legend boxoff
            end
        end
    end

    %% Save all figures
    saveFigures_AC(outdir2)

end



