

function DEV = ERP_DEVPos_AC(events)

ALLEvents = events;

ToneN = 11;
NIter = 1;
NSeq = 1;

while NSeq < 96 && NIter < 7
    ToneS = num2str(ToneN);
    FirstIndi = find(strcmp(ALLEvents, ToneS)); % beginning of the sequence
    NSeq = length(FirstIndi);
    if NSeq < 96
        ToneN = ToneN+10;
    end
    NIter = NIter+1;
end

% Determine position of DEVs

DEV.POS1 = [];
DEV.POS2 = [];
DEV.POS3 = [];
DEV.POS4 = [];

for tt = 1:length(FirstIndi)-1
    
    if find(strcmp(ALLEvents(FirstIndi(tt):FirstIndi(tt+1)-1), '82'))
        DEV.POS1 =  [DEV.POS1 tt];  
    elseif find(strcmp(ALLEvents(FirstIndi(tt):FirstIndi(tt+1)-1), '92'))
        DEV.POS2 =  [DEV.POS2 tt];
    elseif find(strcmp(ALLEvents(FirstIndi(tt):FirstIndi(tt+1)-1), '102'))
        DEV.POS3 =  [DEV.POS3 tt];
    elseif find(strcmp(ALLEvents(FirstIndi(tt):FirstIndi(tt+1)-1), '112'))
        DEV.POS4 =  [DEV.POS4 tt];
    end   
    
end

% And last sequence

if find(strcmp(ALLEvents(FirstIndi(end):end), '82'))
    DEV.POS1 =  [DEV.POS1 tt];
elseif find(strcmp(ALLEvents(FirstIndi(end):end), '92'))
    DEV.POS2 =  [DEV.POS2 tt];
elseif find(strcmp(ALLEvents(FirstIndi(end):end), '102'))
    DEV.POS3 =  [DEV.POS3 tt];
elseif find(strcmp(ALLEvents(FirstIndi(end):end), '112'))
    DEV.POS4 =  [DEV.POS4 tt];
end
    