
function Pack_ERPdata(sess, BSLcorr)


if nargin < 1
    sess = input('Session number?');
    BSLcorr= input('Use BSL-corrected data? 0= No, 1= Yes');
end

% Set up folders

datadir = getSessDir(sess);

if BSLcorr
    outdir1= fullfile(datadir.output, 'ERP', 'BSLcorr');
    files1= dir(fullfile(outdir1, 'blcorr *.mat'));
else
    outdir1= fullfile(datadir.output, 'ERP', 'meancorr'); 
    files1= dir(fullfile(outdir1, 'meancorr *.mat'));
end

load(fullfile(datadir.main, 'ERPheader.mat'));

ResultsDir = fullfile(outdir1, 'Results');
saveName1= fullfile(ResultsDir, 'FullData_MAT.mat');
saveName2= fullfile(ResultsDir, 'DummyData_MAT.mat');

NSubj= datadir.NSubj;
NSeq = datadir.NSeq;

%% Set up 

% Time 

Samp = header.xstep; 
TOI = [-0.1 0.35]; % Time of Interest
Time = header.xstart : Samp : header.xend;
Tind = find(Time >= TOI(1) & Time <= TOI(2));
Time = Time(Tind); % new time vector

% Events of interest

EventCode{1} = {'11', '21', '31', '41', '51', '61', '71', '81', '91', '101'}; % STD
EventCode{2} = {'51', '71', '91'}; % STD strong
EventCode{3} = {'61', '81', '101'}; % STD weak
EventCode{4} = {'82', '92','102','112'}; % DEV
EventCode{5} = {'92','112'}; % DEV s
EventCode{6} = {'82','102'}; % DEV w
EventCode{7} = {'13', '23', '33', '43', '53', '63', '73', '83'}; % STD after 1st DEV

Events = [EventCode{1}, EventCode{4}, EventCode{7}];
EventTypes= length(Events);

% Create GAs for Main Events - This is for the dummy struct

MainEvents = {'STD', 'STDs', 'STDw', 'DEV', 'DEVs', 'DEVw', 'STDreduced'};

Indi{1}(1,:) = find(ismember(Events, EventCode{1}(1:8))); % STD
Indi{2}(1,:) = find(ismember(Events, EventCode{2})); % STDs
Indi{2}(2,:) = find(ismember(Events, EventCode{3})); % STDw
Indi{3}(1,:) = find(ismember(Events, EventCode{4})); % DEV
Indi{4}(1,:) = find(ismember(Events, EventCode{5})); % DEVs
Indi{4}(2,:) = find(ismember(Events, EventCode{6})); % DEVw


%% Here we go

if ~exist(saveName2)
    if ~exist(saveName1)
        
        disp('Creating a new Matrix for Main Events')
        
        %% Load-in and concatenate datasets
       
        disp('Setting up SubjNames, EventCodes, FileNames..')
        
        % Get Subj Names
        
        SubjCode= datadir.SubjNames;
        
        % String name
        
        StringDummy = files1(1).name;
        StringDummy1 = strsplit(StringDummy, '1'); StringDummy1 = StringDummy1{1}; % everything before 'ep'
        StringDummy2 = '.mat';
        
        % Pre-allocate Data mat
        Data= nan(NSubj, EventTypes, NSeq, length(Tind)); % Here we will concatenate Data
        
        for ee= 1:EventTypes
            
            fprintf('\nLooping over Event Code%d\n', ee)
            
            for ss = 1:NSubj
                
                fprintf('\nLooping over Subject%d\n', ss)
                
                % Get file name
                FileName = [StringDummy1, Events{ee}, ' ', SubjCode{ss}, StringDummy2];
                try
                % Load in data
                load(fullfile(outdir1, FileName)); % load specified dataset
                
                % Ch of interest
                
                ChLabels= {'FC chan'};% not in use - but reminder
                chOI= size(data,2); % Select channel of interest
                
                Data(ss,ee,1:size(data,1),:) = squeeze(data(:,chOI,:,:,:,Tind)); % concatenate across subjects
                end
            end
            
            clear data
            
        end
        
        save(saveName1, 'Data', 'Time', 'Events')
        fprintf('\nOk, data are concatenated - moving to dummy structure\n')
        
    else
        load(saveName1)
    end
    
    % Prepare dummy structure with Main Events
    
    dummy = [];
    dummy.time = Time;
    dummy.STD = squeeze(nanmean(nanmean(Data(:,Indi{1}(1,:),:,:),2),3));
    dummy.STDs = squeeze(nanmean(nanmean(Data(:,Indi{2}(1,:),:,:),2),3));
    dummy.STDw = squeeze(nanmean(nanmean(Data(:,Indi{2}(2,:),:,:),2),3));
    dummy.DEV = squeeze(nanmean(nanmean(Data(:,Indi{3}(1,:),:,:),2),3));
    dummy.DEVs = squeeze(nanmean(nanmean(Data(:,Indi{4}(1,:),:,:),2),3));
    dummy.DEVw = squeeze(nanmean(nanmean(Data(:,Indi{4}(2,:),:,:),2),3));
    
    dummy.STDreduced = squeeze(nanmean(nanmean(Data(:,Indi{1}(1,randperm(length(Indi{1}),length(Indi{3}))),:,:),2),3));
    % This has the same N of trials as DEV
    
    save(saveName2, 'dummy')
    
end

