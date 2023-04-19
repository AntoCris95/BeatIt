
function plotTripleBeat_AC(saveName)

% Load in data as provided in inputs

load(saveName); % load(loadName)
PosOI = 1:8;

% Load general beat model 
loadName = strsplit(saveName, 'Triple'); 

% Load Peak mat
load(fullfile(loadName{1}, 'TFR Peaks_MainEvents.mat'))
Peak = Peak.ERS.lowbeta.STD;

% Visualize Triplet models

Names = fieldnames(Models);
NSubj = length(Names);

Labels = {'S-w-w', 'w-S-w', 'w-w-S'};
LNames = {'Triplet1', 'Triplet2', 'Triplet3'};

for ss = 1:NSubj

    % Prepare data of interest from Peak

    data = [];
    data = vertcat(squeeze(Peak(ss*2-1,PosOI,:))', squeeze(Peak(ss*2,PosOI,:))');

    % Display text

    fprintf('\n%s showed a %s pattern \non %d sequences\n', Names{ss}, Labels{1}, length(Triplet1.(Names{ss}).seqN))
    fprintf('\na %s pattern \non %d sequences\n', Labels{2}, length(Triplet2.(Names{ss}).seqN))
    fprintf('\nand a %s pattern \non %d sequences\n', Labels{3}, length(Triplet3.(Names{ss}).seqN))

    Info2Display(ss,1) = round(single(length(Triplet1.(Names{ss}).seqN)),1);
    Info2Display(ss,2) = round(single(length(Triplet2.(Names{ss}).seqN)),1);
    Info2Display(ss,3) = round(single(length(Triplet3.(Names{ss}).seqN)),1);

    % Loop over beat structures

    for l = 1:length(Labels)

        BeatStr = eval(LNames{l});

        % Plot single trials and average

        indi = BeatStr.(Names{ss}).seqN;

        %             fig = figure; fullScreenFig
        %             fig.Name = Names{ss};
        %             subplot(1,3,l); hold on
        %             plot(data(indi,:))

        avgdata = [];
        avgdata = data(indi,:); % keep like this and will average later across subj

        % Prepare GA across subjects per beat-model type
        if ss == 1
            GA.(LNames{l}) = avgdata;
        else
            GA.(LNames{l}) = [GA.(LNames{l}); avgdata];
        end

        % Plot
        %             plot(nanmean(avgdata), 'k','linewidth', 2)
        %
        %             a = gca; Ypos = 0.9*a.YLim(2); Xpos = PosOI(2);
        %             text(Xpos, Ypos, sprintf('%d trials', length(BeatStr.(Names{ss}).seqN)), 'fontweight', 'bold', 'fontsize', 12)
        %             a.XAxis.TickValues = PosOI;
        %             a.FontSize = 12;
        %             a.FontWeight = 'bold';
        %
        %             title(Labels{l})
        %
        %             if l == 1
        %                 xlabel('\fontsize{12} Position along the sequence', 'fontweight', 'bold')
        %                 ylabel('\fontsize{12} Norm amplitude', 'fontweight', 'bold')
        %             end
    end
end
% Save GA
save(saveName, 'Models', 'Triplet1', 'Triplet2', 'Triplet3', 'GA', 'Info2Display')


%% Plot GAs

figure; fullScreenFig

for l = 1:length(Labels)
    
    subplot(1,3,l); hold on
    plot(nanmean(GA.(LNames{l})), 'k','linewidth', 2)
    
    a = gca; Ypos = 0.9*a.YLim(2); Xpos = size(GA.(LNames{1}),2);
    text(Xpos, Ypos, sprintf('%d trials', size(GA.(LNames{l}),2)), 'fontweight', 'bold', 'fontsize', 12)
    a.XAxis.TickValues = PosOI;
    a.FontSize = 12;
    a.FontWeight = 'bold';
    ylabel('% trials')
    
    title(['GA ' Labels{l}])
    
    if l == 1
        xlabel('\fontsize{12} Position along the sequence', 'fontweight', 'bold')
        ylabel('\fontsize{12} Norm amplitude', 'fontweight', 'bold')
    end
end

%% More GA figures and stats

% Plot modelling info

NSeq = 96*2; % sessions
Info2Display = (Info2Display / NSeq) *100;

figure('color', 'white');
Colors = {'b', 'c', 'g', 'm'};
distributionPlot(Info2Display, 'color', Colors); box off
set(gca, 'fontsize', 14, 'fontweight', 'bold', 'xcolor', [1 1 1])
legend(Labels, 'location', 'northwestoutside'); legend boxoff
set(gcf, 'name', 'TripleBeat')

% save Table
writetable(array2table(Info2Display, 'VariableNames', Labels), fullfile(loadName{1}, 'TripleBeat.xls'))

% save figures
path2save = 'C:\Users\anton\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Figures';
saveFigures_AC(fullfile(path2save, 'TripleBeat'))


