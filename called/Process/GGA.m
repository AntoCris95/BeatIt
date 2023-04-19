

%% Create Grand Grand Averages of main conditions from GA data

% This performs Grand-Grand-Average for both ERP and TFR data and
% delivers the Main Events of interest: STD ? DEV * Strong ? Weak.


% Created in March 2020
% Written by Antonio Criscuolo


function GGA(sess, ERPTFR, datatype, BSLcorr)

% Set up

if nargin < 1
    sess = input('Session number?');
    ERPTFR= input('ERP (1) or TFR (2) analysis?');
    datatype= input('Are these simple (0) or complex data (1)?');
    BSLcorr= input('Are these BSL corrected data (1), mean corrected (2) or nocorr (0)?');
end


datadir = getSessDir(sess);

if ERPTFR == 2
    
    if BSLcorr == 1
        outdir= fullfile(datadir.output, 'TFR', 'BSLcorr', 'Results');
    elseif BSLcorr == 2
        outdir= fullfile(datadir.output, 'TFR', 'meancorr', 'Results');
    else
        outdir= fullfile(datadir.output, 'TFR', 'nocorr', 'Results');
    end
    
    lwdata= load(fullfile(datadir.main, 'TFRheader.mat'));
else
    
    if BSLcorr == 1
        outdir= fullfile(datadir.output, 'ERP', 'BSLcorr', 'Results');
    elseif BSLcorr == 2
        outdir= fullfile(datadir.output, 'ERP', 'meancorr', 'Results');
    else
        outdir= fullfile(datadir.output, 'ERP', 'nocorr', 'Results');
    end
    
    lwdata= load(fullfile(datadir.main, 'ERPheader.mat'));
end

if datatype
    files = dir(fullfile(outdir, 'complex GA*.mat'));
    Pos= 3;
else
    files = dir(fullfile(outdir, 'GA *.mat'));
    Pos= 2;
end


MainCond= {'STD', 'STDs', 'STDw', 'DEV', 'DEVs', 'DEVw', ...
    'STDAfterDEV', 'STDsAfterDEV', 'STDwAfterDEV', ...
    '2ndDEV', '2ndDEVs', '2ndDEVw', ...
    'STDafter2DEV'};

%% Check channels

try
    load(fullfile(datadir.main, 'chLabels.mat'));
    chan_labs = chLabels; clear chLabels
catch
    chLabels = 'E:\DATA\Humans - Beatit\DATAproc\chs_OI.mat';
    load(chLabels); %chs_OI
    chan_labs = chs_OI; clear chs_OI
end


%Rremove unwanted chs
Ch2Rem = {'EOGV', 'EOGH', 'A1', 'A2'};
indch = find(ismember(chan_labs,Ch2Rem) < 1);
chan_labs = chan_labs(indch);

Ch= length(chan_labs);

% Fix channels
header= lwdata.header;
for cc = 1:Ch
    header.chanlocs(cc).labels = chan_labs{cc};
end
header.chanlocs(Ch+1:end) = [];


%% Get EventCodes and import Data

fprintf('\nPreparing data for GGA')

Ndatasets= size(files,1);

Data= []; % Here we import Data
EventCode= []; % Here will import list and order of event codes
EventName= []; % Helps us taking out string name

for l= 1:Ndatasets
    EventName= strsplit(files(l).name, ' '); EventCodes{l}= EventName{Pos};
    
    codex= string(EventCodes{l}); codex= strsplit(codex, '.');
    EventCode{l}= codex(1);
    
    load(fullfile(outdir, files(l).name));
    Data= vertcat(Data,data);
end


%% Here we go

for i= 1:length(MainCond)
    
    fprintf('\nFile %d / %d', i, length(MainCond))
    
    % Simulate header structure
        
    if ERPTFR== 2 && datatype
        header.name= ['complex GGA ' MainCond{i}];
    else
        header.name= ['GGA ' MainCond{i}];
    end
    
    MatName= [header.name '.mat'];
    lwName= [header.name '.lw6'];
    
    if ~exist(fullfile(outdir, MatName))
        
        indi= [];
        
        if i == 1 % STD
            
            %         GAEvents= {'31', '41', '51', '61', '71', '81', ...
            %             '91', '101', '111', '121', '131'};
            
            GAEvents= {'51', '61', '71', '81', ...
                '91', '101'};
            
        elseif i == 2 % STDs
            
            %         GAEvents= {'31', '51', '71', '91', '111', '131'};
            
            GAEvents= {'51', '71', '91'};
            
        elseif i == 3 % STDw
            
            %         GAEvents= {'41', '61', '81', '101','121'}; %'121' is noisy
            
            GAEvents= {'61', '81', '101'}; %'121' is noisy
            
        elseif i == 4 % DEV
            
            GAEvents= {'82', '92','102','112'};
            
        elseif i == 5 % DEVs
            
            GAEvents= {'92','112'};
            
        elseif i == 6  % DEVw
            
            GAEvents= {'82','102'};
            
        elseif i == 7 %STD after DEV
            
            GAEvents= {'13', '23', '33', '43', '53', '63', '73', '83'};
            
        elseif i == 8 %STDs after DEV
            
            GAEvents= {'13', '33', '53', '73'};
            
        elseif i == 9 %STDw after DEV
            
            GAEvents= {'23', '43', '63', '83'};
            
        elseif i == 10 %2nd DEV
            
            GAEvents= {'14', '24', '34', '44'};
            
        elseif i == 11 %2nd DEVs
            
            GAEvents= {'14', '34'};
            
        elseif i == 12 %2nd DEVw
            
            GAEvents= {'24', '44'};
            
        elseif i == 13 %STD after 2nd DEV
            
            GAEvents= {'15', '25', '35', '45'};
            
        end
        
        % Go with GGA
        
        for l= 1:length(GAEvents)
            [ind]= find(strcmp(string(EventCode), GAEvents(l)));
            indi(l)= ind;
        end
        
        data= nanmean(Data(indi,:,:,:,:,:),1);
        
        
        if ERPTFR == 2, Time= linspace(lwdata.header.xstart, abs(lwdata.header.xstart), size(Data,6));
        else, Time= linspace(lwdata.header.xstart, lwdata.header.xend, size(Data,6));
        end
        
        header.datasize= size(data);
        Ch = size(data,2);
        header.xstep= round(Time(end) - Time(end-1),3);
        
        header.chanlocs(Ch+1:end)= [];
        header.chanlocs(Ch).labels = 'FC chan';
        
        % Let's get rid of unnecessary events
        
        if ERPTFR == 2
            
            EventName= strsplit(header.name, ' '); EventName= EventName{2};
            
            header.events= [];
            header.events.code= EventName;
            header.events.latency= 0;
            header.events.epoch= 1;
            
        end
        
        % Save avg data
        
        save(fullfile(outdir, MatName), 'data');
        save(fullfile(outdir, lwName), 'header');
        
    end
end

fprintf('\nDONE!\n')
