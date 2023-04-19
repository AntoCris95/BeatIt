
function plotBinary_AC(data, PosOI)

load('C:\Users\p70068941\OneDrive\Work\UM\BAND\Projects\GeneralFunctions\TFR\mycolormap2.mat')

% Calculate amp difference across sequence

data(isnan(data)) = 0;

Ratio = diff(data(:,1:end),[],2); 

% Plot

figure('color', 'white', 'position', [15 63 1878 923]);

subplot(2,1,1); 
plot(Ratio'); box off;
hold on; title('Binary amp diff')
a = gca; a.XAxis.TickValues = PosOI(1:end-1);
a.FontSize = 14;
a.FontWeight = 'bold';

% Calculate pair-wise difference

for trl = 1:size(data,1)
    Similarity(trl,:,:) = bsxfun(@minus,data(trl,:),data(trl,:)');
end

SimMat = squeeze(nanmean(Similarity));
for tl = 1:size(SimMat,1)
    SimMat(tl,tl:end) = nan;
end

% Plot similarity matrix

subplot(2,3,4)
imagesc(abs(SimMat)); colormap(mycolormap2); 
colorbar; box off
a = gca; a.XAxis.TickValues = PosOI;
a.FontSize = 12;
a.FontWeight = 'bold';
title('Pair-wise amp diff')


% Plot grouped similarity

PosInd = [1,3,5];

Mean1 = squeeze(nanmean(nanmean(abs(Similarity(:,PosInd,PosInd)),2),3));
Mean2 = squeeze(nanmean(nanmean(abs(Similarity(:,PosInd+1,PosInd+1)),2),3));

BinarySimilarity = [Mean1;Mean2];

subplot(2,3,5)
distributionPlot([Mean1,Mean2], 'color', {'b', 'c'})
a = gca; 
a.FontSize = 12;
a.FontWeight = 'bold';
a.YLim = [-10 40];

% Go with stats

nperm = 1000;

[p, observeddifference, effectsize, p2Adjust] = permutationTest(Mean1,Mean2, nperm);

% Adjust pvals
[h, p_crit, adj_ci_cvrg, adj_p]= fdr_bh_pCorr(p2Adjust,.05,'pdep');

% Plot as text on figure
if p < p_crit, text2use = 'pAdj < .05';
else, text2use = 'pAdj > .05';  
end   

text(1,a.YLim(2), text2use, 'FontSize', 12, 'FontWeight', 'bold')
title('Binary Similarity')


% Plot dissimilarity

PosInd = [1,3,5];

Mean1 = squeeze(nanmean(nanmean(abs(Similarity(:,PosInd,PosInd+1)),2),3));
Mean2 = squeeze(nanmean(nanmean(abs(Similarity(:,PosInd+1,PosInd)),2),3));

BinaryDissimilarity = [Mean1; Mean2];

subplot(2,3,6)
distributionPlot([BinarySimilarity,BinaryDissimilarity], 'color', {'b', 'c'})
a = gca; 
a.FontSize = 12;
a.FontWeight = 'bold';
a.YLim = [-10 40];

% Go with stats
[p, observeddifference, effectsize, p2Adjust] = permutationTest(BinarySimilarity,BinaryDissimilarity, nperm);
% Adjust pvals
[h, p_crit, adj_ci_cvrg, adj_p]= fdr_bh_pCorr(p2Adjust,.05,'pdep');

% Plot as text on figure
if nanmean(adj_p) < p_crit, text2use = 'pAdj < .05';
else, text2use = 'pAdj > .05';  
end   

text(1,a.YLim(2), text2use, 'FontSize', 12, 'FontWeight', 'bold')
title('Binary Dissimilarity')



%%%%%%%%%%%
% 
% if p < .001, text2use = 'p < .001';
% elseif p == .001, text2use = 'p = .001';    
% elseif p > .05, text2use = 'p > .05';  
% else, text2use = sprintf('pval = %d',round(p,3));
% end   
