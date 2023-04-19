
function modelTripleBeat_AC(saveName, PosOI)

load(saveName); % load(saveName)

BeatInfo = TripleBeat; 

data = GA.TripleBeat;
Subjs = fieldnames(Models); clear Models
NSubj = length(Subjs);

% Prepare for Step-wise regression modelling
% Focus on the triple beat

% 1. Create Predictors: 3 accentuation patters
% Predictor 1: S-w-w

Triplet1 = ones(size(PosOI)) *-.5;
Triplet1(1:3:end) = 1;

% Predictor 2: w-S-w

Triplet2 = ones(size(PosOI)) *-.5;
Triplet2(2:3:end) = 1;

% Predictor 3: w-w-S

Triplet3 = ones(size(PosOI)) *-.5;
Triplet3(3:3:end) = 1;

% 2. Combine predictors

X = [Triplet1; Triplet2; Triplet3];

% 3. Stepwise regression modelling

n = 1; % keeps track of SeqN while retrieving data from the GA
for ss = 1:NSubj

    Name = Subjs{ss}; % SubjN

    % Data to model

    SeqN = BeatInfo.(Name).N; % N of classified seqs
    indi = n:n+SeqN-1; % to retrieve from the GA data

    for tt = 1:SeqN

        y = data(indi(tt),:);

        if ~isnan(y)

            % Ready to model

            Models.(Name).mdl{tt} = stepwiselm(X',y, 'upper', 'linear', 'Criterion', 'adjrsquared', 'verbose', 2); % no interactions

            % Extract pvals
            pVal = Models.(Name).mdl{tt}.coefTest; % pvalue
            Models.(Name).pvals(tt) = pVal;

            % Extract predictors
            preds = Models.(Name).mdl{tt}.PredictorNames;
            Models.(Name).predictors{tt} = preds;

        end
    end

    n = SeqN + 1; % increase index
end

clear Triplet1 Triplet2 Triplet3

% 4. Organize Models and data

for ss = 1:NSubj

t1 = 1; t2 = 1; t3 = 1;

    for seq = 1:length(Models.(Subjs{ss}).predictors)

        if ~isempty(Models.(Subjs{ss}).predictors{seq})
            if strcmp(Models.(Subjs{ss}).predictors{seq}, 'x1') % Pattern 1

                Triplet1.(Subjs{ss}).seqN(t1) = BeatInfo.(Subjs{ss}).seqN(seq);
                t1 = t1+1;

            elseif strcmp(Models.(Subjs{ss}).predictors{seq}, 'x2') % Pattern 2

                Triplet2.(Subjs{ss}).seqN(t2) = BeatInfo.(Subjs{ss}).seqN(seq);
                t2 = t2+1;

            elseif strcmp(Models.(Subjs{ss}).predictors{seq}, 'x3') % Pattern 3

                Triplet3.(Subjs{ss}).seqN(t3) = BeatInfo.(Subjs{ss}).seqN(seq);
                t3 = t3+1;

            end
        end
    end
end

% 5. Save

NewsaveName = strsplit(saveName, 'lowbeta');
NewsaveName = fullfile(NewsaveName{1}, 'TripleBeat_Model.mat');

save(NewsaveName, 'Models', 'Triplet1', 'Triplet2', 'Triplet3')

% 6. Plotting

plotTripleBeat_AC(NewsaveName)
