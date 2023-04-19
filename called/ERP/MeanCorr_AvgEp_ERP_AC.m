

% % Average epochs & apply bsl correction % %

% This function loads in preprocessed and epoched ERP data
% and applies either bsl correction or subtracts mean amplitude or
% expresses data as normalized difference from mean amplitude according to
% what specified.
% Then, it averages epochs. The output is named after the original condition and session.


% Created in March 2020
% Written by Antonio Criscuolo


function MeanCorr_AvgEp_ERP_AC(sess, BSLcorr)


if nargin < 1
    sess= input('Session number?');
    BSLcorr= input('Apply BSL correction? (1), mean correction (2) or no correction (0)?');
end

datadir = getSessDir(sess);
outdir= fullfile(datadir.output, 'ERP');

files = dir(fullfile(outdir, 'ds *.mat')); %Adjust according to the name given to AR files
lwdata= load(fullfile(datadir.main, 'ERPheader.mat'));
    

if BSLcorr == 1
    outdir2 = fullfile(outdir, 'BSLcorr');
elseif BSLcorr == 2
    outdir2 = fullfile(outdir, 'meancorr');
else
    outdir2 = fullfile(outdir, 'nocorr');
end

if ~exist(outdir2)
    mkdir(outdir2);
end


%% Go with data

header= lwdata.header;
load(fullfile(datadir.output, 'chLabels.mat'));
chan_labs = Labels; clear Labels

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

% Channels of interest
chLabels = fullfile(datadir.output, 'chs_OI.mat');
load(chLabels); %chs_OI


for ll= 1:size(files,1)
    
    fprintf('\nFile %d / %d', ll, size(files,1))
    
    load(fullfile(outdir, files(ll).name));
    
    DATA = data; clear data
    
    % Set up ch labels and ROI
    
    Labels = DATA.label;
    Ch= length(Labels);
    
    indch = [];
    for cc = 1:length(chs_OI)
        inds = find(strcmp(Labels, chs_OI{cc}));
        indch = [indch inds];
    end
    
    % Check Events
    
    TotNEvents = length(DATA.events);
    Events = unique(DATA.events);
    NEvents = length(Events);
    
    % Loop over the diff conditions
    
    for ee = 1:NEvents
        
        % Set up file names and folders
        
        name = [Events{ee} ' ' datadir.SubjNames{ll}];
        
        if BSLcorr == 1
            header.name= ['blcorr ' name];
        elseif BSLcorr == 2
            header.name= ['meancorr ' name];
        else
            header.name= ['nocorr ' name];
        end
        
        MatName= [header.name '.mat'];
        lwName= [header.name '.lw6'];
        
        MatName2= ['avg ' header.name '.mat'];
        lwName2= ['avg ' header.name '.lw6'];
        
        if BSLcorr == 1
            saveName= fullfile(outdir, 'BSLcorr', MatName);
            saveHeader= fullfile(outdir, 'BSLcorr', lwName);
            
            saveName2= fullfile(outdir, 'BSLcorr', MatName2);
            saveHeader2= fullfile(outdir, 'BSLcorr', lwName2);
        else
            saveName= fullfile(outdir, 'meancorr', MatName);
            saveHeader= fullfile(outdir, 'meancorr', lwName);
            
            saveName2= fullfile(outdir, 'meancorr', MatName2);
            saveHeader2= fullfile(outdir, 'meancorr', lwName2);
        end
        
        %% Here we go
        
        if ~exist(saveName2)
            if ~exist(saveName)
                
                % Find indices to Event
                ind = find(ismember(DATA.events, Events{ee}));
                
                % Bring to letswave 6D structure
                data = [];
                Fq = 1;

                % Use indices to retrieve data
                for tr = 1:length(ind)
                    data(tr,1:Ch,1,1,1:Fq,:)= DATA.trial{tr};
                end
                
                % Time vector
                time = DATA.time{1};
                
                % Apply Mean / Baseline correction
                
                bsl_corr= ones(size(data)); % pre-allocate
                    
                if BSLcorr == 1
                    
                    % Define bsl
                    % This for ERP
                    bsl= time <= 0 & time >= -0.2;
                    
                    % Apply bsl correction
                    bsl_trial= nanmean(nanmean(data(:,:,:,:,:,bsl),6)); % bsl across all trials
                    bsl_corr= data - bsl_trial;
                    
                elseif BSLcorr == 2
                    
                    % Demean data
                    bsl= time <= 0.3 & time >= -0.2;
                    
                    MeanPow= nanmean(nanmean(data(:,:,:,:,:,bsl),6));
                    bsl_corr= data - MeanPow;
                    
                end
                    
                data = bsl_corr; clear bsl_corr
                data(:,end,:,:,:,:) = nanmean(data(:,indch,:,:,:,:),2);
                
                
                % Adjust header
                
                header.datasize= size(data);
                header.chanlocs(Ch+1:end)= [];
                header.chanlocs(Ch).labels = 'FC chan';
                header.xstart = time(1);
                header.xend = time(end);
                header.xstep = 1/DATA.fsample;
                header.source = datadir.sess;
                
                EventCode= Events{ee};
                
                header.events= [];
                header.events.code= EventCode;
                header.events.latency= 0;
                header.events.epoch= 1;
                
                % Save corrected data
                
                save(saveName, 'data');
                save(saveHeader, 'header');
                
                save(fullfile(datadir.main, 'ERPheader.mat'), 'header');
                
            else
                load(saveName);
            end
            
            %% Calculate mean per event
            
            if size(data,1) > 1
                Data(1,1:Ch,1,1,1,:)= nanmean(data, 1);
                data = Data;
            end

            header.datasize= size(data);
            header.name = ['avg ' header.name];
            
            % Save avg data
            
            save(saveName2, 'data');
            save(saveHeader2, 'header');
            
            clear data Data
            
        end
    end
end
