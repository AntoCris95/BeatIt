


function Pack_TFR_FullData(sess, BSLcorr)


% Set up

if nargin < 1
    sess = input('Session number?');
    BSLcorr= input('Use BSL-corrected data? 0= No, 1= Yes');
end

datadir = getSessDir(sess);
datadir.output = fullfile(datadir.sess, 'DATAproc');

if BSLcorr == 1
    outdir= fullfile(datadir.output, 'TFR', 'BSLcorr');
    files = dir(fullfile(outdir, 'bl*.mat'));
else
    outdir= fullfile(datadir.output, 'TFR', 'meancorr');
    files = dir(fullfile(outdir, 'meancorr*.mat'));
end

lwdata= load(fullfile(datadir.main, 'TFRheader.mat'));
header = lwdata.header; clear lwdata

ResultsDir = fullfile(outdir, 'Results');
saveName= fullfile(ResultsDir, 'FullData_MAT.mat');

NSubj= datadir.NSubj;
NSeq = datadir.NSeq;

% Time

Samp = 1/250;
TOI = [-0.4 0.5]; % Time of Interest
Time = -1 : Samp : 1; % Time - after segmentation
TimeOrig = -2 : Samp : 2; % Time - before segmentation
Tind = find(Time >= TOI(1) & Time <= TOI(2));
Time = Time(Tind); % new time vector

% Define fq of interest

FQS = [1 40];
FQStep= round(1 / (abs(TimeOrig(1)) + TimeOrig(end)),3); % Frequency resolution = inverse of the period
Fqs = FQS(1) : FQStep : FQS(2);

% Events of interest

MainEvents = {'STD', 'STDs', 'STDw', 'DEV', 'DEVs', 'DEVw'};

EventCode{1} = {'11', '21', '31', '41', '51', '61', '71', '81', '91', '101'}; % STD
EventCode{2} = {'51', '71', '91'}; % STD strong
EventCode{3} = {'61', '81', '101'}; % STD weak
EventCode{4} = {'82', '92','102','112'}; % DEV
EventCode{5} = {'92','112'}; % DEV s
EventCode{6} = {'82','102'}; % DEV w
EventCode{7} = {'13', '23', '33', '43', '53', '63', '73', '83'}; % STD after 1st DEV

Events = [EventCode{1}, EventCode{4}];


%% Here we go

if ~exist(saveName)
    
    disp('Creating a new Matrix for Main Events')
    
    %% Load-in and concatenate datasets
    
    disp('Setting up SubjNames, EventCodes, FileNames..')
    
    % Get Subj Names
    
    SubjCode= datadir.SubjNames;
    
    % String name
    
    StringDummy = files(1).name;
    StringDummy1 = strsplit(StringDummy, ' '); StringDummy1 = StringDummy1{1}; % 'meancorr'
    StringDummy2 = '.mat';
    
    % Load and pack data according to eventsOI
    
    for i = 1:length(MainEvents) % loop over main events
        
        EventTypes= length(EventCode{i});
        
        % Pre-allocate Data mat
        Data = [];
        Data = nan(NSubj, EventTypes, NSeq, length(Fqs), length(Tind)); % Here we will concatenate Data
        
        for ee= 1:EventTypes
            
            fprintf('\nLooping over Event Code%d\n', ee)
            
            for ss = 1:NSubj
                
                fprintf('\nLooping over Subject%d\n', ss)
                
                % Get file name
                FileName = [StringDummy1, ' ', Events{ee}, ' ', SubjCode{ss}, StringDummy2];
                
                % Load in data
                load(fullfile(outdir, FileName)); % load specified dataset
                
                % Ch of interest
                ChLabels= {'FC chan'};% not in use - but reminder
                chOI= size(data,2); % Select channel of interest
                
                Data(ss,ee,1:size(data,1),:,:) = single(squeeze(data(:,chOI,:,:,:,Tind))); % concatenate across subjects
                
            end
            clear data
        end
        
        data = single(Data); 
        clear Data
        
        NElements = size(data,2);
        
        % Let's re-organize data according to frequencies:
        
        %Delta
        
        delta1= [];
        delta1.time= Time;
        Fq= find(Fqs >= 1.5 & Fqs <= 4);
        delta1.nanmean= squeeze(nanmean(data(:, :, :, Fq, :),4));
        delta1.SE= (std(squeeze(nanmean(data(:, :, :, Fq, :),4)), 'omitnan'))/sqrt(length(Fq)+NElements);
        
        %Theta
        
        theta1= [];
        theta1.time= Time;
        Fq= find(Fqs >= 4 & Fqs <= 8);
        theta1.nanmean= squeeze(nanmean(data(:, :, :, Fq, :),4));
        theta1.SE= (std(squeeze(nanmean(data(:, :, :, Fq, :),4)), 'omitnan'))/sqrt(length(Fq)+NElements);
        
        %Alpha
        
        alpha1= [];
        alpha1.time= Time;
        Fq= find(Fqs >= 9 & Fqs <= 13);
        alpha1.nanmean= squeeze(nanmean(data(:, :, :, Fq, :),4));
        alpha1.SE= (std(squeeze(nanmean(data(:, :, :, Fq, :),4)), 'omitnan'))/sqrt(length(Fq)+NElements);
        
        %Beta
        
        beta1= [];
        beta1.time= Time;
        Fq= find(Fqs >= 12 & Fqs <= 25);
        beta1.nanmean= squeeze(nanmean(data(:, :, :, Fq, :),4));
        beta1.SE= (std(squeeze(nanmean(data(:, :, :, Fq, :),4)), 'omitnan'))/sqrt(length(Fq)+NElements);
        
        %Low-Beta
        
        lowbeta1= [];
        lowbeta1.time= Time;
        Fq= find(Fqs >= 12 & Fqs <= 20);
        lowbeta1.nanmean= squeeze(nanmean(data(:, :, :, Fq, :),4));
        lowbeta1.SE= (std(squeeze(nanmean(data(:, :, :, Fq, :),4)), 'omitnan'))/sqrt(length(Fq)+NElements);
        
        %High-Beta
        
        highbeta1= [];
        highbeta1.time= Time;
        Fq= find(Fqs >= 20 & Fqs <= 25);
        highbeta1.nanmean= squeeze(nanmean(data(:, :, :, Fq, :),4));
        highbeta1.SE= (std(squeeze(nanmean(data(:, :, :, Fq, :),4)), 'omitnan'))/sqrt(length(Fq)+NElements);
        
        
        if i == 1
            
            STD = [];
            STD.delta= delta1;
            STD.theta= theta1;
            STD.alpha= alpha1;
            STD.beta= beta1;
            STD.lowbeta= lowbeta1;
            STD.highbeta= highbeta1;
            
        elseif i == 2
            
            STDs = [];
            STDs.delta= delta1;
            STDs.theta= theta1;
            STDs.alpha= alpha1;
            STDs.beta= beta1;
            STDs.lowbeta= lowbeta1;
            STDs.highbeta= highbeta1;
            
        elseif i == 3
            
            STDw = [];
            STDw.delta= delta1;
            STDw.theta= theta1;
            STDw.alpha= alpha1;
            STDw.beta= beta1;
            STDw.lowbeta= lowbeta1;
            STDw.highbeta= highbeta1;
            
        elseif i == 4
            
            DEV = [];
            DEV.delta= delta1;
            DEV.theta= theta1;
            DEV.alpha= alpha1;
            DEV.beta= beta1;
            DEV.lowbeta= lowbeta1;
            DEV.highbeta= highbeta1;
            
        elseif i == 5
            
            DEVs = [];
            DEVs.delta= delta1;
            DEVs.theta= theta1;
            DEVs.alpha= alpha1;
            DEVs.beta= beta1;
            DEVs.lowbeta= lowbeta1;
            DEVs.highbeta= highbeta1;
            
        elseif i == 6
            
            DEVw = [];
            DEVw.delta= delta1;
            DEVw.theta= theta1;
            DEVw.alpha= alpha1;
            DEVw.beta= beta1;
            DEVw.lowbeta= lowbeta1;
            DEVw.highbeta= highbeta1;
            
        end
    end
    
    % Save
    save(saveName, MainEvents{:});
    
    
else
    disp('Data structure ready for plots')
end
