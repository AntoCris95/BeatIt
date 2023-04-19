
function go_FastTFR_AC(sess)

%% Set up

if nargin < 1
    sess= input('Session number?');
end

datadir = getSessDir(sess);
outdir= fullfile(datadir.sess, 'DATAproc','TFR');

load(fullfile(datadir.main, 'TFRheader.mat'));
files = dir(fullfile(outdir, 'ds *.mat')); 

Fq_OI= [1 40];


%% Here we go

load(fullfile(datadir.main, 'chLabels.mat'));
chan_labs = Labels; clear Labels

%Rremove unwanted chs
Ch2Rem = {'EOGV', 'EOGH', 'A1', 'A2'}; 
indch = find(ismember(chan_labs,Ch2Rem) < 1);
chan_labs = chan_labs(indch);

Ch= length(chan_labs);

% Fix channels
for cc = 1:Ch
    header.chanlocs(cc).labels = chan_labs{cc};
end
header.chanlocs(Ch+1:end) = [];

% Channels of interest
chLabels = 'E:\DATA\Humans - Beatit\DATAproc\chs_OI.mat';
load(chLabels); %chs_OI


for ff = 1:length(files)
    
    fprintf('\n Preparing file %d / %d for TFR', ff, length(files))
    
    load(fullfile(outdir, files(ff).name));
    
    Time = data.time{1};
    FQStep= round(1 / (abs(Time(1)) + Time(end)),3); % Frequency resolution = inverse of the period
    FqAll = Fq_OI(1) : FQStep : Fq_OI(2);
    
    % Check Events
    
    ALLEvents = data.events;
    TotNEvents = length(data.events);
    Events = unique(data.events);
    NEvents = length(Events);
    
    DATA = data;
    
    fprintf('\n Using parallel computing for TFR')
            
    for ee = 1:NEvents
        
        fprintf('\n Looping over Event %d / %d', ee, NEvents)
        
        trials2use = find(ismember(ALLEvents, Events{ee}));
        saveName= ['fftcwt ' Events{ee} ' ' files(ff).name];
        
        if ~exist(fullfile(outdir, saveName))
            
            Data = zeros([length(trials2use), size(DATA.trial{1},1), size(DATA.trial{1},2)]); % Trials x Chs X Time
            
            for tt = 1:length(trials2use)
                for cc = 1:size(DATA.trial{1},1)
                    Data(tt,cc,:) = DATA.trial{trials2use(tt)}(cc,:);
                end
            end
            
            
            % % Initiate parallel loops % %
            try, parpool, end
            
            % Iterate through channels %
            
            TF_Data = zeros([size(Data,1), size(Data,2), length(FqAll), size(Data,3)]);
            
            parfor cc = 1:size(Data,2)
                TF_Data(:,cc,:,:) = FFTWavelet_AC(squeeze(Data(:,cc,:)), Time, Fq_OI)
            end
            
            
            % Prepare for letswave data structure
            
            % Create FC chan
            
            Labels = DATA.label;
            
            indch = [];
            
            for cc = 1:length(chs_OI)
                inds = find(strcmp(Labels, chs_OI{cc}));
                indch = [indch inds];
            end
            
            clear data Data % letswave will need 'data' when saving
            
            % Reduce time dimension
            
            TOF = [-1 1];
            tind = find(Time >= TOF(1) & Time <= TOF(2));
            
            Data = [];
            Data = zeros([size(TF_Data,1), size(TF_Data,2)+1, size(TF_Data,3), length(tind)]);
            Data = TF_Data(:,:,:,tind);
            Data(:,end+1,:,:) = nanmean(TF_Data(:,indch,:,tind),2);
            
            
            % Bring into letswave format
            
            Nrep= size(Data, 1);
            Nch= size(Data, 2);
            Nfqs= size(Data, 3);
            Ntpoints= size(Data, 4);
            
            data(1:Nrep,1:Nch,1,1,1:Nfqs,1:Ntpoints)= Data;
            data = single(data);
            
            % Simulate header structure for letswave
            
            headName= strsplit(saveName, '.');  headName= headName{1};
            header.name = headName;
            header.datasize = size(data);
            header.xstart = TOF(1);
            header.ystart = Fq_OI(1);
            header.xstep = Time(2) - Time(1);
            header.xend = TOF(2);
            
            header.chanlocs(Nch+1:end)= [];
            header.chanlocs(Nch).labels = 'FC chan';
            
            header.events= [];
            header.events.code= Events{ee};
            header.events.latency= 0;
            header.events.epoch= 1;
            
            %% Save
            
            fprintf('\n And saving!')
            
            save(fullfile(outdir, saveName), 'data', '-v7.3');
            save(fullfile(outdir, [headName '.lw6']), 'header');
            
            clear data Data
            
        end
    end
end



