

% This function runs preprocessing using a combination of
% custom-made functions and scripts and the Fieldtrip toolbox.

% In this preprocessing pipeline, we perform, in order:
% 1. Import raw data in *.cnt/.hdr format;
% 2. Preprocessing: band-pass frequency filtering, rereference to mastoids;
% 3. ICA, by using fastica algorithm. This part has a semiautomatic (visual
% inspection) and fully automatic options;
% 4. Segmentation: import original .*trg files, edit trigger info according
% to specific hypothesis and use it to create epochs. 
% Epoch length  = -2 to 2 time-locked to stimulus onset
% 5. Additional round of artifact correction:
% Two options: in a semi-automatic procedure, we go into visual inspection of data,
% and you can specificy whether there's need for further artifact corrections.
% Thus, two more options:
% A. standard artifact rejection (AR): trial rejection based on amplitude criteria;
% B. artifact suppression (AS): artifacts are identified based on a thresh criterion and then
% we substitue noisy time-windows by means of temporal cubic interpolation.
% In the fully automatic version, we directly go into AR/AS.
% Such step ensures (almost) no trials are lost;
% 6. Call ft_preprocessing again: 
% A. Second round of frequency-filtering: parameters differ for ERP and TFR data;
% B. Linear detrending and remove high-order trends (polyremoval; both optional)
% NB: ERP data are re-segmented to -.5 to .5s relative to tone onsets
% before the new ft_preprocessing.
% 7. Down-sampling: half the original sampling frequency.


% Data folders and other information are retrieved by calling getSessDir
% function. Specify sessN, whether you want to work with ERP or TFR data
% and which pipeline you want to use.

% Created in March 2021
% Written by Antonio Criscuolo


function Preprocess_AC(sess, ERPTFR, STDpipeline)

if nargin < 2
    sess = input('Session number?');
    ERPTFR = input('ERP (1) or TFR (2) analysis?');
    STDpipeline = input('Pipeline with artifact rejection (1) or suppression(2)?'); % 2= artifact suppression with cubic interpolation
end

datadir = getSessDir(sess);
outdir= fullfile(datadir.sess, 'DATAproc');
Normalize = 1; % 0 = no normalization, works with a threshold of 85uV; 1 = normalize, threshold between 3-5. 

if ERPTFR == 1
    outdir2 = fullfile(outdir, 'ERP');
else
    outdir2 = fullfile(outdir, 'TFR');
end

if ~exist(outdir)
    mkdir(outdir)
end
if ~exist(outdir2)
    mkdir(outdir2)
end

files = dir(fullfile(datadir.sess, '*.hdr'));


% Here we go

AutoICA = 1; % automatic ICA = 1; manual = 0;

% Loop over files

for ff = 1:length(files)
    
    % Set up names first
    
    rawfilename= fullfile(datadir.sess, files(ff).name);
    SubjName = strsplit(files(ff).name, '.'); SubjName = SubjName{1};
    SubjName = [SubjName '.mat'];
    
    % Step 1: import and preprocess
    
    saveName= fullfile(outdir, SubjName);
    
    Name2save1 = ['but1 reref ' SubjName];
    saveName1= fullfile(outdir, Name2save1);
    
    % Step 2: apply independent-component analysis
    
    Name2save2 = ['ICA ' Name2save1];
    saveName2= fullfile(outdir, Name2save2);
    
    % Step 3: segmentation
    
    Name2save3 = ['ep ' Name2save2];
    saveName3= fullfile(outdir, Name2save3);
    
    % Step 4: artifact correction/suppression
    
     if STDpipeline == 1
        Name2save4 = ['AR ' Name2save3];
    else
        Name2save4 = ['AS ' Name2save3];
        Name2save4andAR = ['AR ' Name2save4];
     end
    
    saveName4= fullfile(outdir, Name2save4);
    saveName4andAR= fullfile(outdir, Name2save4andAR);
    
    % Step 5: new band-pass filtering and detrending
    
    Name2save5 = ['but2 nodtrend ' Name2save4]; %temporary 'nodetrend'
    saveName5= fullfile(outdir2, Name2save5);
    
    % Step 6: downsample
    
    Name2save6 = ['ds ' Name2save5];
    saveName6= fullfile(outdir2, Name2save6);
    
    
    % Ready to go: let's check whether files already exist
    
    if ~exist(saveName6) % ds
        if ~exist(saveName5) % new band-pass fq filter and no detrending
            if ~exist(saveName4) % AR or AS
                if ~exist(saveName3) % segmentation & trend removal
                    if ~exist(saveName2) % ICA
                        if ~exist(saveName1) % initial preproc
                            if ~exist(saveName) % import into mat
                                
                                % Import raw data
                                
                                fprintf('\nImporting dataset N%d\n', ff)
                                
                                data = read_refa8_data(rawfilename);
                                save(saveName, 'data');
                                
                            else
                                load(saveName);
                            end
                            
                            % Get Labels
                            
                            if ff == 1
                                for ll = 1:length(data.label)-1
                                    Labels{ll} = data.label{ll};
                                end
                                save(fullfile(datadir.sess, 'chLabels.mat'), 'Labels');
                            else
                                load(fullfile(datadir.sess, 'chLabels.mat'));
                            end
                            
                            % Preprocess data
                            
                            fprintf('\nPreprocessing dataset N%d\n', ff)
                            
                            cfg=[];
                            cfg.channel = Labels;
                            cfg.demean = 'no';
                            cfg.continuous = 'yes';
                            cfg.dftfilter = 'no';
                            cfg.bpfilter = 'yes'; % band-pass; other options: lp/hp/bs filters
                            cfg.bpfreq = [0.1 50]; % according to what's up: bp/lp/hp freq
                            cfg.bpfiltord = 4;
                            cfg.detren = 'no'; % yes
                            cfg.reref = 'yes';
                            cfg.refchannel = {'A1', 'A2'};
                            
                            % Distributed computing %
                            cfg.inputfile= saveName;
                            cfg.outputfile= saveName1; % saves automatically
                            
                            data = ft_preprocessing(cfg);
                            
                        else
                            load(saveName1);
                        end
                        
                        
                        % % Perform ICA and remove components % %
                        % Will use distributed computing within FieldTrip
                        % ICA data is automatically saved
                        
                        fprintf('\nApplying ICA on dataset N%d\n', ff)
                        
                        data = goICA_AC(data, outdir, Name2save1);
                        
                    else
                        load(saveName2)
                    end
                    
                    % Segmentation
                    
                    fprintf('\nSegmentation on dataset N%d\n', ff)
                    
                    data = Segmentation_AC(data, datadir, ff);
                    events = data.events;
                    
                    save(saveName3, 'data');
                    
                else
                    load(saveName3);
                    events = data.events;
                end
                
                
                % Artifact rejection / suppression
                
                selchan = ft_channelselection({'all', '-A1', '-A2'}, data.label);
                cfg = [];
                cfg.channel = selchan;
                data = ft_selectdata(cfg, data);
                
                if STDpipeline == 1
                    
                    fprintf('\nPerforming artifact rejection on dataset N%d\n', ff)
                    
                    [data, Nartifacts] = ArtRejection_AC(data, ERPTFR);
                    
                else
                    
                    if AutoICA == 0
                        
                        % Check how it looks like
                        
                        cfg = [];
                        cfg.viewmode = 'vertical';
                        cfg.continuous = 'no';
                        cfg.artifactalpha = 0.8;
                        cfg.blocksize = 5;
                        ft_databrowser(cfg, data);
                        
                        NeedAS= input('Need further Artifact suppression? 1 (yes), 0 (no)');
                        
                    else
                        NeedAS = 1;
                    end
                    
                    close all
                    
                    if NeedAS == 1
                        
                        fprintf('\nPerforming artifact suppression on dataset N%d\n', ff)
                        
                        % Perform artifact suppression by means of temporal linear interpolation.
                        % Computations are run ch-by-ch and trial-by-trial for accuracy.
                        % Will make use of parallel computing toolbox to speed it up.
                        
                        [data, Artifacts] = ArtSuppression_AC(data, Normalize);
                        
                        save(saveName4, 'data');
                        
                        % Double-check: remove trials with unsuccessful AS
                        [data, Nartifacts2] = ArtRejection_AC(data);
                        events = data.events;
                        
                    end
                    
                    SaveName = ['Nart ' datadir.SubjNames{ff} '.mat'];
                    savingName= fullfile(outdir, SaveName);
                    
                    save(saveName4andAR, 'data');
                    try % there may be no additional artifacts
                        save(savingName, 'Artifacts', 'Nartifacts2');
                    catch
                        save(savingName, 'Artifacts');
                    end
                end
                
            else
                load(saveName4andAR);
                events = data.events;
            end
            
            % Shorten ERP data

            if ERPTFR == 1
                TOI = [-.5 .5];
                data = redefineTrial_AC(data, TOI);
            end

            % New fq filtering
            
            fprintf('\nNew round of Fq filtering on dataset N%d\n', ff)
            
            cfg = [];
            cfg.channel = 'EEG';
            cfg.demean = 'no'; % whole trial
            cfg.continuous = 'no';
            cfg.dftfilter = 'no';
            if ERPTFR == 1
                cfg.bpfilter = 'yes'; % band-pass; other options: lp/hp/bs filters
                cfg.bpfreq = [1 30]; % according to what's up: bp/lp/hp freq
                cfg.bpfiltord = 4;
            else
                cfg.lpfilter = 'yes'; % low-pass; other options: lp/hp/bp filters
                cfg.lpfreq = 40; % according to what's up: bp/lp/hp freq
                cfg.lpfiltord = 4;
            end
            cfg.detren = 'no';
            cfg.reref = 'no';
            cfg.polyremoval = 'no';
            
            data = ft_preprocessing(cfg, data);
            data.events = events;
            
            save(saveName5, 'data');
            
        else
            load(saveName5);
            events = data.events;
        end
        
        % Downsampling
        
        fprintf('\nResampling dataset N%d\n', ff)
        
        cfg = [];
        cfg.resamplefs = data.fsample/2;
        data = ft_resampledata(cfg, data);
        data.events = events;
        
        save(saveName6, 'data');
        
    end
end


