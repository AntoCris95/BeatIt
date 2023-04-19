
function TFR_Peak2Peak_Sequence2(sess, BSLcorr)

%% Peak-to-Peak Amplitude and Latency

% This function loads in preprocessed (either bsl- or not-bsl-corrected;
% and either ERP of TFR) data, and detects peak amplitude and latecy for
% every event of interest in three time-windows corresponding to
% event-related P50, N100 and P200.
% Then, it generates GAs for all Events.
% These info are stored in a Peak matrix which can be later used to generate
% plots and stats.


% Edited in April 2021
% Written by Antonio Criscuolo


if nargin < 1
    sess = input('Session number?');
    BSLcorr= input('Use BSL-corrected data? 0= No, 1= Yes');
end


% Set up folders

datadir = getSessDir(sess);

if BSLcorr
    outdir1= fullfile(datadir.output, 'TFR', 'BSLcorr');
    files1= dir(fullfile(outdir1, 'blcorr *.mat'));
else
    outdir1= fullfile(datadir.output, 'TFR', 'meancorr');
    files1= dir(fullfile(outdir1, 'meancorr *.mat'));
end

ResultsDir = fullfile(outdir1, 'Results');
lwdata= load(fullfile(datadir.main, 'TFRheader.mat'));
loadName = fullfile(ResultsDir, 'FullData_MAT.mat');
saveName= fullfile(ResultsDir, 'TFR Peaks_MainEvents.mat');

NSubj= datadir.NSubj;

%% Load data mat

if ~exist(loadName)
    Pack_TFR_FullData(sess, BSLcorr)
end

%% Find peak, amplitude and latency for time-points of interest

% We focus on two event-related components:
% 1: ERD. Draw a window of +- 100ms
% 2: ERS. Draw a window of +- 100ms

% On these frequency-bands:
% 1. Delta, 1-4Hz;
% 2. Theta, 4-8Hz;
% 3. Alpha, 8-12Hz
% 4. Lowbeta, 12-20Hz
% 5. Beta, 12-25Hz
% 6. Highbeta, 20-25Hz

FqOI = {'delta', 'theta', 'alpha', 'lowbeta', 'beta', 'highbeta'};
CentFq = [2, 6, 10, 16, 19, 22]; % center frequencies

% And 2 main conditions:
% 1. STD
% 2. DEV


if ~exist(saveName)
    
    disp('Creating a new Peak Matrix for Main Events')
    
    Peak= []; % This will store the matrix
    load(loadName)
    
    fprintf('\nDefining time-windows and proceding to Peaks extraction\n')
    
    Time = STD.alpha.time;
    
    % Loop over frequencies
    % Use dynamic structure fields
    
    for fqs = 1:length(FqOI)
        
        % Define Time frame of interest
        TimePad = round(1/CentFq(fqs),2); % symmetrical padding proportional to center frequency
        if TimePad > .15 % but not too much
            TimePad = .15;
        end
        
        % ERD
        % Find peak
        Event= -.15;
        TOF= [Event-TimePad 0-TimePad]; 
        TimeIndERD= find(Time >= TOF(1) & Time <= TOF(2));
        
        [~,tindi] = min(nanmean(STD.(FqOI{fqs}).nanmean(:,:,:,TimeIndERD)), [],4);
        
        tindi = squeeze(tindi);
%         TimePad = TimePad/2;
            
        % Redefine time
        
        for pos = 1:size(tindi,1)
            for seq = 1:size(tindi,2)
                
                Event= Time(TimeIndERD(tindi(pos,seq)));
                TOF= [Event-TimePad Event+TimePad];
                TimeIndERD2= find(Time >= TOF(1) & Time <= TOF(2));
                
                % Create Peak structure for ERD
                
                [Peak.ERD.(FqOI{fqs}).STD(:,pos,seq)]= nanmean(STD.(FqOI{fqs}).nanmean(:,pos,seq,TimeIndERD2),4);
                if pos <= size(DEV.(FqOI{fqs}).nanmean,2)
                [Peak.ERD.(FqOI{fqs}).DEV(:,pos,seq)]= nanmean(DEV.(FqOI{fqs}).nanmean(:,pos,seq,TimeIndERD2),4);
                end
            end
        end
        
        % ERS
        % Find peak
        
%         TimePad = TimePad*2;
        
        Event= .1;
        TOF= [0+TimePad Event+TimePad];
        TimeIndERS= find(Time >= TOF(1) & Time <= TOF(2));
        
        [~,tindi] = max(nanmean(STD.(FqOI{fqs}).nanmean(:,:,:,TimeIndERS)),[],4);
        
        tindi = squeeze(tindi);
%         TimePad = TimePad/2;
        
        % Redefine time
        
        for pos = 1:size(tindi,1)
            for seq = 1:size(tindi,2)
                
                Event= Time(TimeIndERS(tindi(pos,seq)));
                TOF= [Event-TimePad Event+TimePad];
                TimeIndERS2= find(Time >= TOF(1) & Time <= TOF(2));
                
                % Create Peak structure for ERS
                
                [Peak.ERS.(FqOI{fqs}).STD(:,pos,seq)]= nanmean(STD.(FqOI{fqs}).nanmean(:,pos,seq,TimeIndERS2),4);
                if pos <= size(DEV.(FqOI{fqs}).nanmean,2)
                [Peak.ERS.(FqOI{fqs}).DEV(:,pos,seq)]= nanmean(DEV.(FqOI{fqs}).nanmean(:,pos,seq,TimeIndERS2),4);
                end
            end
        end
    end
    
    

    %% Save Peak structure
    
    fprintf('\nAnd saving!\n')
    
    save(saveName, 'Peak');
    
else
    
    disp('Peak Matrix already exist');
    load(saveName)
    
end

%% Plot

% Plot_Peak2Peak_Sequence(Peak)

