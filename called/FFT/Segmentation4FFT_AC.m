
% This function is called during the newPreprocess_AC pipeline
% and applies segmentation to data according to specified parameters.

% Use data, datadir, name2save, fileN and ERPTFR as input parameters.
% Trigger information are automatically retrieved from raw *.trg files,
% edited according to experimental hypothesis and used to determine
% conditions and trials.
% Sample info are then used in the following artifact correction procedure.

% Created in March 2021
% Written by Antonio Criscuolo


function data = Segmentation4FFT_AC(data, datadir, fileN)

% Time info

Samp = 1/data.fsample;
Tpoints = length(data.time{1});
Time = data.time{1};

% Setup parameters

Ntones = 13; % min 13 tones; NB: we will exclude first 1s 
ISO = 0.65; %ms
StimRate = round(1/0.65, 2); %Hz

Start = 0; % Start from the beginning of the sequence
% Start = 2*ISO; % start from the third tone (NB: the first is at 0)
PostStim = round(Ntones/StimRate); % will be around 8s

TpointsStart = abs(Start) / Samp;
TpointsEnd = PostStim / Samp;

SeqTime = (1:Ntones)*ISO - ISO; % From first onset (0) til 13th tone onset

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
    DEVpos(tr) = cell2mat( textscan(alltext{2}, '%d'));
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
% Apply segmentation

if ~isempty(EventsOI)
    
    % B. Extract onsets
    for tr = 1:length(EventsOI)
        Onsets(tr) = triggers(EventsOI(tr)).sample;
    end
    
    NEvents = length(EventsOI); % N of sequences
    
    % Every trigger code signals the beginning of a sound sequence
    % Define segments 
    
%     Onset = Onsets - TpointsPre; % if starting before seq
    Onset = Onsets + TpointsStart; % if starting bit after the seq
    Offset = Onset + TpointsEnd;
    
    % Create trials
    
    for ee = 1:NEvents % N of sequences
%         try % we may miss the last sound sequence
            
            SegData.trial{ee} = data.trial{1}(:,Onset(ee) : Offset(ee));
            SegData.time{ee} = linspace(Start, PostStim, length(SegData.trial{ee}));
            SegData.sampleinfo(ee,:) = [Onset(ee), Offset(ee)];
            
            % Add event labels
            
            % STD
            STDlabels = {'STD'}; STDlabels = repmat(STDlabels, [1,DEVpos(ee)-1]);
            SegData.events{ee} = STDlabels;
            
            % DEV
            SegData.events{ee}{end+1} = 'DEV';
            
            % STD after DEV - sometimes 2 DEV
            STDlabels = {'STDpostDEV'};
            
            if DEVnum(ee) == 2
                DEVposi = 12;
            else
                DEVposi = DEVpos(ee);
            end
            Npos = SEQnum(ee) - DEVposi;
            SegData.events{ee}(DEVposi+1:SEQnum(ee)) = repmat(STDlabels, [1,Npos]);
%         end
    end
    
else
    fprintf('Stopped at fileN %d', fileN)
    error('Could not find Events of interest')
end

data.trial = SegData.trial;
data.time = SegData.time;
data.sampleinfo = SegData.sampleinfo;
data.Ntones = Ntones;
data.seqStart = Start;
data.seqEnd = PostStim;
data.stimRate_Hz = StimRate;
data.ISO = ISO;
data.events = SegData.events;



