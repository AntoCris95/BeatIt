
%% Create GA for event type

% This is the Grand-Average function for both ERP and TFR data.
% Output is named 'GA %s' per event code.

% Created in March 2020
% Written by Antonio Criscuolo


function GA(sess, ERPTFR, datatype, BSLcorr)

%Set up

if nargin < 1
    sess = input('Session number?');
    ERPTFR= input('ERP (1) or TFR (2) analysis?');
    if ERPTFR == 2
        datatype= input('Are these simple (0) or complex data (1)?');
    end
    BSLcorr= input('Are these BSL corrected data (1), mean corrected (2) or nocorr (0)?');
end


datadir = getSessDir(sess);
datadir.output = fullfile(datadir.sess, 'DATAproc');

if ERPTFR == 2
    
    if BSLcorr == 1
        outdir= fullfile(datadir.output, 'TFR', 'BSLcorr');       
    else
        outdir= fullfile(datadir.output, 'TFR', 'meancorr');
    end
            
    if datatype == 1
        files = dir(fullfile(outdir, 'avg nocorr cwtcomplex*.mat'));
    else
        files = dir(fullfile(outdir, 'avg *.mat'));
    end
    
    lwdata= load(fullfile(datadir.main, 'TFRheader.mat'));
    
else
    
    if BSLcorr == 2
        outdir= fullfile(datadir.output, 'ERP', 'meancorr');
        files = dir(fullfile(outdir, 'avg meancorr *.mat'));
    elseif BSLcorr == 1
        outdir= fullfile(datadir.output, 'ERP', 'BSLcorr');
        files = dir(fullfile(outdir, 'avg blcorr *.mat'));
    else
        outdir= fullfile(datadir.output, 'ERP', 'nocorr');
        files = dir(fullfile(outdir, 'avg nocorr *.mat'));
    end
    
    lwdata= load(fullfile(datadir.main, 'ERPheader.mat'));
    
end

if ~exist(fullfile(outdir, 'Results'))
    mkdir(fullfile(outdir, 'Results'))
end

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

%% Here we go

disp('Preparing data for GA')

Ndatasets= size(files,1);

% Define position of eventName in the header.name string

EventName= strsplit(files(1).name, ' ');
Pos= length(EventName)-1;

% Retrieve event codes

for ff = 1:Ndatasets
    event = strsplit(files(ff).name, ' ');
    EventCode{ff} = event{Pos};
end

% Determine event types

ff = 1; ee = 1;
while ff < Ndatasets
    Events{ee}= EventCode{ff};
    [ind] = find(strcmp(EventCode{ff}, EventCode));
    ff = ind(end)+1; ee = ee+1;
end

EventTypes= length(Events);

% Concatenate datasets

for i=1:EventTypes
    
    fprintf('\nFile %d / %d', i, EventTypes)
    
    % Simulate header structure
    
    header= lwdata.header;
    
    if ERPTFR == 2 && datatype
        header.name= ['complex GA ' Events{i}];
    else
        header.name= ['GA ' Events{i}];
    end
    
    MatName= [header.name '.mat'];
    lwName= [header.name '.lw6'];
    
    if ~exist(fullfile(outdir, 'Results', MatName))
        
        Data= [];
        
        [ind]= find(strcmp(EventCode, Events{i}));
        
        for l=1:length(ind)
            load(fullfile(outdir, files(ind(l)).name));
            Data= vertcat(Data,nanmean(data,1));
        end
        
        clear data
        
        % Calculate mean per event code
        
        data= nanmean(Data,1);
        
        if ERPTFR == 2, Time= linspace(lwdata.header.xstart, abs(lwdata.header.xstart), size(Data,6));
        else, Time= linspace(lwdata.header.xstart, lwdata.header.xend, size(Data,6));
        end
        
        header.datasize= size(data);
        Ch = size(data,2);
        header.xstep= round(Time(end) - Time(end-1),3);
        
        header.chanlocs(Ch+1:end)= [];
        header.chanlocs(Ch).labels = 'FC chan';
        
        
        % Let's get rid of unnecessary events
        
        header.events= [];
        header.events.code= Events{i};
        header.events.latency= 0;
        header.events.epoch= 1;
        
        
        % Save avg data
        
        save(fullfile(outdir, 'Results', MatName), 'data');
        save(fullfile(outdir, 'Results', lwName), 'header');
        
    end
end

