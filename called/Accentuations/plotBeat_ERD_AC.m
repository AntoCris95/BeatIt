
function plotBeat_ERD_AC(data2load1, data2load2, FqOI, PosOI)

%% Load in data as provided in inputs

load(data2load1); % load(loadName)
load(data2load2) % load(saveName)
saveName = data2load2;
NPos = length(PosOI);

%% Visualize models

Names = fieldnames(Models);
NSubj = length(Names);

Labels = {'BinaryBeat', 'TripleBeat', 'CombBeat', 'OtherBeat'};

if ~exist('GAERD', 'var')
    for ss = 1:NSubj
        
        % Display text
        fprintf('\n%s showed a binary beat pattern \non %d sequences\n', Names{ss}, BinaryBeat.(Names{ss}).N)
        fprintf('\na triplet beat pattern \non %d sequences\n', TripleBeat.(Names{ss}).N)
        fprintf('\nand a combined binary + triplet beat pattern \non %d sequences\n', CombBeat.(Names{ss}).N)
        fprintf('\n%d sequences were not classified\n', OtherBeat.(Names{ss}).N)
        
        Info2Display(ss,1) = BinaryBeat.(Names{ss}).N;
        Info2Display(ss,2) = TripleBeat.(Names{ss}).N;
        Info2Display(ss,3) = CombBeat.(Names{ss}).N;
        Info2Display(ss,4) = OtherBeat.(Names{ss}).N;
        
        
        % Data to retrieve
        y = Peak.ERD.(FqOI).STD(ss:ss+1,PosOI,:);
        
        % Organize data
        data = [];
        for sess = 1:2
            tempdata = squeeze(y(sess,:,:));
            tempdata = reshape(tempdata, [numel(tempdata) 1]);
            data = vertcat(data, tempdata);
        end
        
        % Prepare figure
        fig = figure('color', 'white', 'position', [16 238 1878 625]);
        fig.Name = Names{ss};
        
        % Loop over beat structures
        
        for l = 1:length(Labels)
            
            BeatStr = eval((Labels{l}));
            
            subplot(1,4,l); hold on
            
            % Create trial avg for plotting
            avgdata = [];
            for nn = 1:BeatStr.(Names{ss}).N
                indi = BeatStr.(Names{ss}).seqN(nn)*NPos-NPos+1 : BeatStr.(Names{ss}).seqN(nn)*NPos;
                plot(data(indi))
                avgdata = [avgdata; data(indi)'];
            end
            
            % Prepare GA across subjects per beat-model type
            if ss == 1
                GAERD.(Labels{l}) = avgdata;
            else
                GAERD.(Labels{l}) = [GAERD.(Labels{l}); avgdata];
            end
            
            % Plot
            plot(nanmean(avgdata), 'k','linewidth', 2)
            
            a = gca; Ypos = 0.9*a.YLim(2); Xpos = PosOI(2);
            text(Xpos, Ypos, sprintf('%d trials', BeatStr.(Names{ss}).N), 'fontweight', 'bold', 'fontsize', 12)
            a.XAxis.TickValues = PosOI;
            a.FontSize = 12;
            a.FontWeight = 'bold';
            
            title(Labels{l})
            
            if l == 1
                xlabel('\fontsize{12} Position along the sequence', 'fontweight', 'bold')
                ylabel('\fontsize{12} Norm amplitude', 'fontweight', 'bold')
            end
        end
    end
    % Save GA
    save(saveName, 'Models', 'BinaryBeat', 'TripleBeat', 'CombBeat', 'OtherBeat', 'GA', 'Info2Display', 'GAERD')
end

%% Plot GAs

fig = figure('color', 'white', 'position', [16 238 1878 625]);

for l = 1:length(Labels)
    
    subplot(1,4,l); hold on
    %     plot(GA.(Labels{l})');
    plot(nanmean(GAERD.(Labels{l})), 'k','linewidth', 2)
    
    a = gca; Ypos = 0.9*a.YLim(2); Xpos = PosOI(2);
    text(Xpos, Ypos, sprintf('%d trials', length(GAERD.(Labels{l}))), 'fontweight', 'bold', 'fontsize', 12)
    a.XAxis.TickValues = PosOI;
    a.FontSize = 12;
    a.FontWeight = 'bold';
    
    title(['GA ERD ' Labels{l}])
    
    if l == 1
        xlabel('\fontsize{12} Position along the sequence', 'fontweight', 'bold')
        ylabel('\fontsize{12} Norm amplitude', 'fontweight', 'bold')
    end
end

%% Extract amplitude ratio across positions

plotBinary_AC(GAERD.BinaryBeat,PosOI)

plotBinary2_AC(GAERD.OtherBeat,PosOI)





