

% % Average epochs & apply bsl correction % %

% This function loads in preprocessed and epoched TFR data
% and applies either bsl correction or subtracts mean amplitude or
% expresses data as normalized difference from mean amplitude according to
% what specified.
% Then, it averages epochs. The output is named after the original condition and session.
% Bsl correction is expressed as %change in case of TFR amplitude;
% in case of complex TFR, there is no correction.

% Created in March 2021
% Written by Antonio Criscuolo


function MeanCorr_AvgEp_TFR_AC(sess, datatype, BSLcorr)


if nargin < 1
    sess = input('Session number?');
    BSLcorr = input('Apply BSL correction? (1), mean correction (2) or no correction (0)?');
    datatype = input('Are these simple (0) or complex data (1)?');
end

datadir = getSessDir(sess);
datadir.output = fullfile(datadir.sess, 'DATAproc');

outdir= fullfile(datadir.output, 'TFR');

if datatype
    files = dir(fullfile(outdir, 'cwtcomplex *.mat'));
    BSLcorr = 0; disp('Not applying correction to complex data')
else
    files = dir(fullfile(outdir, 'fftcwt *.mat')); %or cwt
end

lwdata= load(fullfile(datadir.main, 'TFRheader.mat'));

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

load(fullfile(datadir.main, 'chLabels.mat'));
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


for ll= 1:size(files,1)
    
    fprintf('\nFile %d / %d', ll, size(files,1))
          
    % Check Events
    
    Events= strsplit(files(ll).name, ' ');
    Events= Events{2};
    
    % Set up file names and folders
    
    SubjName = strsplit(files(ll).name, '.');
    SubjName = strsplit(SubjName{1}, ' ');
    SubjName = SubjName{end};
    
    name = [Events ' ' SubjName];
    
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
    
    saveName= fullfile(outdir2, MatName);
    saveHeader= fullfile(outdir2, lwName);

    saveName2= fullfile(outdir2, MatName2);
    saveHeader2= fullfile(outdir2, lwName2);
    
    %% Here we go
    
    if ~exist(saveName2)
        if ~exist(saveName)
            
            load(fullfile(outdir, files(ll).name));
            
            % We're already in letswave 6D format
            
            Fq = size(data,5);
            time = linspace(-1, 1, size(data, 6));
            
            % Apply Mean / Baseline correction
            
            bsl_corr= ones(size(data)); % pre-allocate
            
            if datatype == 0 && BSLcorr == 1
                
                % This for TFR
                bsl= time <= 0 & time >= -0.2;
                
                bsl_trial= nanmean(nanmean(data(:,:,:,:,:,bsl),6));
                bsl_corr= ((data - bsl_trial)./ (data + bsl_trial)) *100;
                
            elseif BSLcorr == 2 && datatype == 0
                
                % Demean data
                bsl= time <= 0.3 & time >= -0.2;
                
                % Express data as normalized difference from mean amplitude
                MeanPow= nanmean(nanmean(data(:,:,:,:,:,bsl),6));
                bsl_corr= ((data - MeanPow)./ (data + MeanPow)) *100;
                MeanPow = nanmean(nanmean(bsl_corr(:,:,:,:,:,bsl),6));
                bsl_corr = bsl_corr - MeanPow;
                
            else
                
                % For complex exponentials, subtraction of the
                % imaginary part would not be intepretable here
                bsl_corr= data;
                
            end
            
            data = bsl_corr; clear bsl_corr
            
            
            % Adjust header
            
            header.datasize= size(data);
            header.chanlocs(Ch+1).labels = 'FC chan';
            header.ystart = 1;
            header.xstart = time(1);
            header.xend = time(end);
            header.xstep = time(end)-time(end-1);
            header.source = datadir.sess;
            
            EventCode= Events;
            
            header.events= [];
            header.events.code= EventCode;
            header.events.latency= 0;
            header.events.epoch= 1;
            
            % Save corrected data
            
            save(saveName, 'data', 'header', '-v7.3');
            save(saveHeader, 'header');
            
            save(fullfile(datadir.main, 'TFRheader.mat'), 'header');
            
        else
            load(saveName);
        end
        
        %% Calculate mean per event
        
        Size = size(data); Size(1) = 1;
        Data = zeros(Size);
        if size(data,1) > 1
            Data= nanmean(data, 1);
        end
        data = Data; clear Data
        
        header.datasize= size(data);
        header.name = ['avg ' header.name];
        
        % Save avg data
        
        save(saveName2, 'data', '-v7.3');
        save(saveHeader2, 'header');
        
        clear data Data
        
    end
    
end
