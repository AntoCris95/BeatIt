
% This function is called during the newPreprocess_AC pipeline
% and applies segmentation to data according to specified parameters.

% Use data, datadir, name2save, fileN and ERPTFR as input parameters.
% Trigger information are automatically retrieved from raw *.trg files,
% edited according to experimental hypothesis and used to determine
% conditions and trials.
% Sample info are then used in the following artifact correction procedure.

% Created in March 2021
% Written by Antonio Criscuolo


function data = Segmentation_AC(data, datadir, fileN, ERPTFR)


% Setup parameters

Samp = 1/data.fsample; % sampling

ISI = 0.6; % inter-stimulus interval: 600ms
Dur = 0.05; % tone duration: 50ms

Tpoints = (ISI+Dur) / Samp;

% Pre - Post stimulus

if nargin == 4
    if ERPTFR == 1
        PreStim = -0.5;
        PostStim = 0.5; % 0.8 for longer seq
    else
        PreStim = -2;
        PostStim = 2;
    end
else
    PreStim = -2;
    PostStim = 2;
end
    

TpointsPre = abs(PreStim) / Samp;
TpointsPost = PostStim / Samp;


% Read-in trig files

% 1. Txt file:
% Stores info about sequence lenght and DEV numbers and position

Trg = dir(fullfile(datadir.sess, '*.txt'));
TrgName = fullfile(datadir.sess, Trg(fileN).name);

trigInfo = readtable(TrgName);
trigCode = table2cell(trigInfo(:,end));
trigInfo = table2cell(trigInfo(:,5));

% Extract info
for tr = 1:size(trigInfo,1)
    alltext = strsplit(char(trigInfo(tr)), '_');
    DEVnum(tr) = cell2mat(textscan(alltext{end}, '%d'));
    DEVpos(tr) = cell2mat(textscan(alltext{2}, '%d'));
    SEQnum(tr) = cell2mat(textscan(alltext{3}, '%d'));
end

%2. Hdr file:
% Stores triggers info

trigfiles = dir(fullfile(datadir.sess, '*.hdr'));
triggers = read_refa8_event(fullfile(datadir.sess, trigfiles(fileN).name));

% Returns Nx1 struct array with N triggers
%   trg(i).type     ... type
%   trg(i).sample   ... linear index
%   trg(i).value    ... trigger code (double)
%   trg(i).time     ... trigger latency in s
%   trg(i).offset   ... byte offset
%   trg(i).duration ... duration

% Edit trig file:
% A. Select conds of interest

for ee = 1:length(triggers)
    TriggersALL(ee) = triggers(ee).value;
end

EventsOI = find(ismember(TriggersALL, cell2mat(trigCode)));


% If trigger codes match between hdr and trig file, go on

if ~isempty(EventsOI)
    
    % B. Extract onsets
    for tr = 1:length(EventsOI)
        Onsets(tr) = triggers(EventsOI(tr)).sample;
    end
    
    % Apply segmentation
    
    NEvents = length(EventsOI);
    
    % Every trigger code signals the beginning of a sound sequence
    % Thus, we need to calculate the length of each sequence based on the
    % N of tones. Next, we label DEV tones.
    % Loop over the trigger codes and further loop over each tone in the
    % sequence.
    
    TrialN = 0;
    
    for ee = 1:NEvents
        
        NTones = SEQnum(ee);
        NDev = DEVnum(ee);
        
        for tt = 1:NTones % sequence order
            
            nt = TrialN + tt; % trialN, takes into account all iterations
            
            TimeStart = Onsets(ee) + Tpoints*(tt-1) - TpointsPre;
            TimeEnd = Onsets(ee) + Tpoints*(tt-1) + TpointsPost;
            
            SegData.sampleinfo(nt,:) = [TimeStart TimeEnd];
            SegData.trial{nt} = data.trial{1}(:,TimeStart : TimeEnd);
            SegData.time{nt} = linspace(PreStim, PostStim, length(SegData.trial{nt}));
            
            if tt < DEVpos(ee)
                
                SegData.event{nt} = [num2str(tt) '1']; % STD before DEV
                
            elseif tt == DEVpos(ee)
                
                SegData.event{nt} = [num2str(tt) '2']; % 1st DEV
                
            elseif tt > DEVpos(ee) && NDev == 1
                
                SegData.event{nt} = [num2str(tt-DEVpos(ee)) '3']; % STD after DEV
                
            elseif tt > DEVpos(ee) && NDev == 2 && tt < 12
                
                SegData.event{nt} = [num2str(tt-DEVpos(ee)) '3']; % STD after DEV
                
            elseif NDev == 2 && tt == 12 
                
                SegData.event{nt} = [num2str(tt-DEVpos(ee)) '4']; % 2nd DEV
               
            elseif NDev == 2 && tt > 12
                     
                SegData.event{nt} = [num2str(tt-12) '5']; % STD after 2nd DEV
                     
            end
        end
        
        TrialN = TrialN + tt; % update num Trials
        
    end
    
    data.trial = SegData.trial;
    data.time = SegData.time;
    data.events = SegData.event;
    data.sampleinfo = double(SegData.sampleinfo);
    
else
    fprintf('Stopped at SubjN %d', fileN)
    error('Trial codes do not match!')
    return
end

