

function modelBeat_AC(sess, BSLcorr)


if nargin < 1
    sess = input('Session number?');
    BSLcorr= input('Use BSL-corrected data? 0= No, 1= Yes');
end

% Set up parameters of interest

% FQS = {'Alpha', 'Lowbeta', 'Beta', 'Highbeta'}; % old peak structure
% FQS = {'delta', 'theta', 'alpha', 'lowbeta', 'beta', 'highbeta'};

FQS = {'lowbeta'}; % single fq call for plotting

PosOI = 1:8; % positions of interest along the auditory sequence
NPos = length(PosOI);

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
loadName= fullfile(ResultsDir, 'TFR Peaks_MainEvents.mat');

NSubj= datadir.NSubj;


%% Here we go

for fqs = 1:length(FQS)
    
    saveName = fullfile(ResultsDir, sprintf('%s RegrModels.mat', FQS{fqs}));
    FqOI = FQS{fqs};
    
    if ~exist(saveName)
        
        % Load data mat with peak amps
        
        if ~exist(loadName)
%             TFR_Peak2Peak_Sequence(sess, BSLcorr)
            TFR_Peak2Peak_Sequence2(sess, BSLcorr)
        end
        load(loadName)
        
        
        % Prepare for Step-wise regression modelling
        
        DivFact = 1; % change amplitude - not in use
        n = 1; % keeps track of the subject number
        
        for ss = 1:2:NSubj % pool together 2sessions per subject
            
            % Data to model
            y = Peak.ERS.(FqOI).STD(ss:ss+1,PosOI,:);
            
            newdata = [];
            
            for sess = 1:2
                tempdata = squeeze(y(sess,:,:));
                tempdata = reshape(tempdata, [numel(tempdata) 1]);
                newdata = vertcat(newdata, tempdata);
            end
            
            y = newdata;
            
            % Predictor 1 = Binary beat
            BinaryModel = ones(size(y));
            BinaryModel(2:2:end) = -1; % S-w pattern
            BinaryModel = BinaryModel / DivFact;
            
            % Predictor 2 = Triplet beat
            TripletModel = ones(size(y))*-0.5;
            TripletModel(1:3:end) = 1; % S-w-w pattern
            %             TripletModel(2:3:end) = 0; % S-w-w pattern
            TripletModel = TripletModel / DivFact;
            
            % Predictor 3 = Error term
            %         ErrorTerm = reshuffle(y .* randn(size(y)));
            ErrorTerm = rand(size(y));
            
            % Predictor 4 = Constant
            Constant = ones(size(y));
            
            % Combine predictors
            
            %              X = [BinaryModel TripletModel Constant ErrorTerm];
            X = [BinaryModel TripletModel Constant];
            
            
            % Go with stepwise regression modelling - loop over trials
            
            Name = sprintf('Subj%d', n);
            
            SeqN = 1; % keeps track of the sequence number
            for tt = 1:NPos:length(y)
                
                indi = tt:tt+NPos-1; % position indices

                if ~isnan(y(indi))
                
                % Model
                Models.(Name).mdl{SeqN} = stepwiselm(X(indi,:),y(indi), 'upper', 'linear', 'Criterion', 'adjrsquared', 'verbose', 2); % no interactions
                
                % Extract pvals
                pVal = Models.(Name).mdl{SeqN}.coefTest; % pvalue
                Models.(Name).pvals(SeqN) = pVal;
                
                % Extract predictors
                preds = Models.(Name).mdl{SeqN}.PredictorNames;
                Models.(Name).predictors{SeqN} = preds;

                end
                SeqN = SeqN + 1; % increase index
            end
            n = n+1; % increase index
        end
        
        
        % Summarize models
        
        Names = fieldnames(Models);
        NSubj = length(Names);
        
        for ss = 1:NSubj
            
            bb = 1; % binary beat
            tb = 1; % triplet beat
            cb = 1; % combined beat
            ob = 1; % other beat
            
            % Find significant models
            %         indi = find(Models.(Names{ss}).pvals <= .05);
            
            % Scan through all models
            indi = 1:length(Models.(Names{ss}).predictors);
            
            for ii = 1:length(indi)
                
                if length(Models.(Names{ss}).predictors{indi(ii)}) % non-empty string
                    if strcmp(Models.(Names{ss}).predictors{indi(ii)}(1), 'x1') % only binary or binary + error term
                        
                        try % combined beat
                            if strcmp(Models.(Names{ss}).predictors{indi(ii)}(2), 'x2') % = binary + triplet
                                CombBeat.(Names{ss}).seqN(cb) = indi(ii);
                                CombBeat.(Names{ss}).N = cb;
                                CombBeat.(Names{ss}).Preds{cb} = Models.(Names{ss}).predictors{indi(ii)};
                                cb = cb+1;
                            else
                                BinaryBeat.(Names{ss}).seqN(bb) = indi(ii);
                                BinaryBeat.(Names{ss}).N = bb;
                                BinaryBeat.(Names{ss}).AdjRsquared(bb) = Models.(Names{ss}).mdl{indi(ii)}.Rsquared.Adjusted;
                                BinaryBeat.(Names{ss}).Rsquared(bb) = Models.(Names{ss}).mdl{indi(ii)}.Rsquared.Adjusted;
                                BinaryBeat.(Names{ss}).Bcoeff(bb) = Models.(Names{ss}).mdl{indi(ii)}.Coefficients(2,1);
                                bb = bb+1;
                            end
                        catch
                            BinaryBeat.(Names{ss}).seqN(bb) = indi(ii);
                            BinaryBeat.(Names{ss}).N = bb;
                            BinaryBeat.(Names{ss}).AdjRsquared(bb) = Models.(Names{ss}).mdl{indi(ii)}.Rsquared.Adjusted;
                            BinaryBeat.(Names{ss}).Rsquared(bb) = Models.(Names{ss}).mdl{indi(ii)}.Rsquared.Adjusted;
                            BinaryBeat.(Names{ss}).Bcoeff(bb) = Models.(Names{ss}).mdl{indi(ii)}.Coefficients.Estimate(2);
                            bb = bb+1;
                        end
                        
                    elseif strcmp(Models.(Names{ss}).predictors{indi(ii)}(1), 'x2') % only triplet or triplet + error term
                        TripleBeat.(Names{ss}).seqN(tb) = indi(ii);
                        TripleBeat.(Names{ss}).N = tb;
                        TripleBeat.(Names{ss}).AdjRsquared(tb) = Models.(Names{ss}).mdl{indi(ii)}.Rsquared.Adjusted;
                        TripleBeat.(Names{ss}).Rsquared(tb) = Models.(Names{ss}).mdl{indi(ii)}.Rsquared.Adjusted;
                        TripleBeat.(Names{ss}).Bcoeff(tb) = Models.(Names{ss}).mdl{indi(ii)}.Coefficients.Estimate(2);
                        tb = tb+1;
                    else
                        OtherBeat.(Names{ss}).seqN(ob) = indi(ii);
                        OtherBeat.(Names{ss}).N = ob;
                        ob = ob+1;
                    end
                else
                    OtherBeat.(Names{ss}).seqN(ob) = indi(ii);
                    OtherBeat.(Names{ss}).N = ob;
                    ob = ob+1;
                end
            end
            
        end
        save(saveName, 'Models', 'BinaryBeat', 'TripleBeat', 'CombBeat', 'OtherBeat')
    end
end


%% Plot Binary beat models on Beta ERS

PosOI = 1:8;

disp('Calling plotting function')
plotBeat_AC(loadName, saveName, FqOI, PosOI)

% Plot modelling info

load(saveName)

figure('color', 'white');
Colors = {'b', 'c', 'g', 'm'};
distributionPlot(Info2Display, 'color', Colors)
set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1])

%% Plot Binary beat models on Beta ERD

plotBeat_ERD_AC(loadName, saveName, FqOI, PosOI)

%% Focus on the Triple beat now

modelTripleBeat_AC(saveName, PosOI) % also calls the plotting function



