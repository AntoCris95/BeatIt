
function ShadedPlots_TFR_CompleteData(sess, BSLcorr)


% Set up

if nargin < 1
    sess = input('Session number?');
    BSLcorr= input('Use BSL-corrected data? 0= No, 1= Yes');
end

datatype= 0; disp('Using simple (0) data ');
datadir = getSessDir(sess);
datadir.output = fullfile(datadir.sess, 'DATAproc');

if BSLcorr == 1
    
    outdir= fullfile(datadir.output, 'TFR', 'BSLcorr');
    
    if datatype == 1
        files = dir(fullfile(outdir, 'bl cwtcomplex*.mat'));
    else
        files = dir(fullfile(outdir, 'bl cwt*.mat'));
    end
    
else
    
    outdir= fullfile(datadir.output, 'TFR', 'meancorr');
    
    if datatype == 1
        files = dir(fullfile(outdir, 'avg nocorr *.mat'));
    else
        files = dir(fullfile(outdir, 'avg meancorr*.mat'));
    end
    
end

lwdata= load(fullfile(datadir.main, 'TFRheader.mat'));
header = lwdata.header; clear lwdata

saveName= fullfile(outdir, 'Results', 'Data_MAT.mat');

% Define data of interest

Fq2Plot = 6; % delta, theta, alpha, beta, high-beta, low-beta

Comparisons = {'STD s-w', 'DEV s-w', 'STD-DEV'};
PlannedContrasts(1,:) = [2,3];
PlannedContrasts(2,:) = [5,6];
PlannedContrasts(3,:) = [1,4];

Labels= {'STD', 'STDs', 'STDw', 'DEV', 'DEVs', 'DEVw'};


% Load in data

Ndatasets= size(files,1);

% Define position of eventName in the header.name string

EventName= strsplit(files(1).name, ' ');
Pos= length(EventName)-1;

% Retrieve event codes

for ff = 1:Ndatasets
    event = strsplit(files(ff).name, ' ');
    EventCode{ff} = event{Pos};
end

% Determine event types

ff = 1; ee = 1;
while ff < Ndatasets
    Events{ee}= EventCode{ff};
    [ind] = find(strcmp(EventCode{ff}, EventCode));
    ff = ind(end)+1; ee = ee+1;
end

EventTypes= length(Events);


% Concatenate datasets of interest

if ~exist(saveName)
    disp('Preparing TFR data for plots')
    Pack_TFRdata(sess, BSLcorr)
    load(saveName)
else
    load(saveName)
end


%% Let's plot

disp('Preapring shaded plots')

% Order

figOrd(1,:)= 1:Fq2Plot; % subplots: alpha, beta, lowbeta, highbeta
figOrd(2,:)= Fq2Plot+1:Fq2Plot*2;
figOrd(3,:)= Fq2Plot*2+1:Fq2Plot*3; % three figures as the N of contrasts

% Single plots

for i = 1:length(Comparisons)
    for ee = 1:2 % two elements per plot
        
        n = 1;
        
        if i == 1 && ee == 1
            
            delta1 = STDs.delta;
            theta1 = STDs.theta;
            alpha1 = STDs.alpha;
            beta1= STDs.beta;
            lowbeta1 = STDs.lowbeta;
            highbeta1 = STDs.highbeta;
            
        elseif i == 1  && ee == 2
            
            delta1 = STDw.delta;
            theta1 = STDw.theta;
            alpha1 = STDw.alpha;
            beta1= STDw.beta;
            lowbeta1 = STDw.lowbeta;
            highbeta1 = STDw.highbeta;
        
        elseif i == 2  && ee == 1
            
            delta1 = DEVs.delta;
            theta1 = DEVs.theta;
            alpha1 = DEVs.alpha;
            beta1= DEVs.beta;
            lowbeta1 = DEVs.lowbeta;
            highbeta1 = DEVs.highbeta;
        
        elseif i == 2  && ee == 2
            
            delta1 = DEVw.delta;
            theta1 = DEVw.theta;
            alpha1 = DEVw.alpha;
            beta1= DEVw.beta;
            lowbeta1 = DEVw.lowbeta;
            highbeta1 = DEVw.highbeta;
            
        elseif i == 3 && ee == 1
            
            delta1 = STD.delta;
            theta1 = STD.theta;
            alpha1 = STD.alpha;
            beta1= STD.beta;
            lowbeta1 = STD.lowbeta;
            highbeta1 = STD.highbeta;
            
        elseif i == 3  && ee == 2
            
            delta1 = DEV.delta;
            theta1 = DEV.theta;
            alpha1 = DEV.alpha;
            beta1= DEV.beta;
            lowbeta1 = DEV.lowbeta;
            highbeta1 = DEV.highbeta;
            
        end

        if ee == 1, Colors = '-b'; else Colors = '-r'; end
        
        figure(figOrd(i,n)); hold on
        box off, set(gcf,'color','w')
        shadedErrorBar(delta1.time, delta1.nanmean, delta1.SE, Colors); %'lineProps', 
        title([Comparisons{i} ' Delta avg ch'])
        
        plot1=gcf; plot1= plot1.CurrentAxes;
        plot1.FontSize = 12;
        plot1.FontWeight = 'bold';
        xlabel('\fontsize{12} Time (s)', 'fontweight', 'bold')
        ylabel('\fontsize{12} Norm % change', 'fontweight', 'bold')
        %     legend(Labels)
        n = n+1;
        
        figure(figOrd(i,n)); hold on
        box off, set(gcf,'color','w')
        shadedErrorBar(theta1.time, theta1.nanmean, theta1.SE, Colors); %'lineProps', 
        title([Comparisons{i} ' Theta avg ch'])
        
        plot1=gcf; plot1= plot1.CurrentAxes;
        plot1.FontSize = 12;
        plot1.FontWeight = 'bold';
        xlabel('\fontsize{12} Time (s)', 'fontweight', 'bold')
        ylabel('\fontsize{12} Norm % change', 'fontweight', 'bold')
        %     legend(Labels)
        n = n+1;
        
        figure(figOrd(i,n)); hold on
        box off, set(gcf,'color','w')
        shadedErrorBar(alpha1.time, alpha1.nanmean, alpha1.SE, Colors); %'lineProps',
        title([Comparisons{i} ' Alpha avg ch'])
        
        plot1=gcf; plot1= plot1.CurrentAxes;
        plot1.FontSize = 12;
        plot1.FontWeight = 'bold';
        xlabel('\fontsize{12} Time (s)', 'fontweight', 'bold')
        ylabel('\fontsize{12} Norm % change', 'fontweight', 'bold')
        %     legend(Labels)
        n = n+1;
        
        figure(figOrd(i,n)); hold on
        box off, set(gcf,'color','w')
        shadedErrorBar(beta1.time, beta1.nanmean, beta1.SE, Colors); %'lineProps', 
        title([Comparisons{i} ' Beta avg ch'])
        
        plot1=gcf; plot1= plot1.CurrentAxes;
        plot1.FontSize = 12;
        plot1.FontWeight = 'bold';
        xlabel('\fontsize{12} Time (s)', 'fontweight', 'bold')
        ylabel('\fontsize{12} Norm % change', 'fontweight', 'bold')
        %     legend(Labels)
        n = n+1;
        
        figure(figOrd(i,n)); hold on
        box off, set(gcf,'color','w')
        shadedErrorBar(lowbeta1.time, lowbeta1.nanmean, lowbeta1.SE, Colors); %'lineProps', 
        title([Comparisons{i} ' Lowbeta avg ch'])
        
        plot1=gcf; plot1= plot1.CurrentAxes;
        plot1.FontSize = 12;
        plot1.FontWeight = 'bold';
        xlabel('\fontsize{14} Time (s)', 'fontweight', 'bold')
        ylabel('\fontsize{14} Norm % change', 'fontweight', 'bold')
        %     legend(Labels)
        n = n+1;
        
        figure(figOrd(i,n)); hold on
        box off, set(gcf,'color','w')
        shadedErrorBar(highbeta1.time, highbeta1.nanmean, highbeta1.SE, Colors); %'lineProps', 
        title([Comparisons{i} ' Highbeta avg ch'])
        
        plot1=gcf; plot1= plot1.CurrentAxes;
        plot1.FontSize = 12;
        plot1.FontWeight = 'bold';
        xlabel('\fontsize{12} Time (s)', 'fontweight', 'bold')
        ylabel('\fontsize{12} Norm % change', 'fontweight', 'bold')
        %     legend(Labels)
        n = n+1;
        

    end
end

% Focus on STD s-w with stats

lowbeta1 = STDs.lowbeta;
lowbeta2 = STDw.lowbeta;

figure; set(gcf, 'position', [267 213 667 430], 'color', 'white'); hold on

h1 = shadedErrorBar(lowbeta1.time, lowbeta1.nanmean, lowbeta1.SE, 'lineProps','b', 'transparent',1); %'lineProps',
h2 = shadedErrorBar(lowbeta1.time, lowbeta2.nanmean, lowbeta2.SE, 'lineProps','r', 'transparent',1); %'lineProps',

title([Comparisons{1} ' Lowbeta avg ch'])
AdjAxes_AC; ylabel('Norm % change')
leg = legend([h1.mainLine, h2.mainLine], Labels{2:3}); legend boxoff


% Statistical testing

TimeOI = [0 .2];
Tind = lowbeta1.time >= TimeOI(1) & lowbeta1.time <= TimeOI(2);

% figure; hold on
% plot(STD.Time, nanmean(lowbeta1.all))
% plot(STD.Time, nanmean(lowbeta2.all))

alpha = .05; %alpha/nperm
[H,P,CI,STATS] = ttest(lowbeta1.all(:,Tind), lowbeta2.all(:,Tind));

% Adjust pvals
[h, p_crit, adj_ci_cvrg, adj_p]= fdr_bh_pCorr(P,.05,'pdep');

% Plot stats
Stars = nan(size(adj_p));
Stars2 = nan(size(P));

Stars(find(adj_p < alpha)) = 1.3;
Stars2(find(P < alpha)) = 1.3;
plot(lowbeta1.time(Tind), Stars, '*k');
plot(lowbeta1.time(Tind), Stars2, '*r');

leg.String{3} = 'pAdj <.05';
leg.String{4} = 'p <.05';

clear p observeddifference effectsize p2Adjust Data1 data
