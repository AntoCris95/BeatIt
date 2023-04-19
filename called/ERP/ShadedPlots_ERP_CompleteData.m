
function ShadedPlots_ERP_CompleteData(sess, BSLcorr)

% Create plots with shaded error bars


% Created in September 2020
% Written by Antonio Criscuolo


% Set up

if nargin < 1
    sess = input('Session number?');
    BSLcorr= input('Use BSL-corrected data? 0= No, 1= Yes');
end

if BSLcorr == 0, disp('Using non-bsl corrected data'),
else disp('Using bsl-corrected data'), end


% Set up folders

datadir = getSessDir(sess);

if BSLcorr == 1
    outdir= fullfile(datadir.output, 'ERP', 'BSLcorr');
else
    outdir= fullfile(datadir.output, 'ERP', 'meancorr');
end
load(fullfile(datadir.main, 'ERPheader.mat'));

Fqs= 1;
NSubj = datadir.NSubj;

dataName= fullfile(outdir, 'Results', 'DummyData_MAT.mat');

% Load in data

disp('Loading in data')
try load(dataName); dummy
catch
    Pack_ERPdata(sess, BSLcorr)
    load(dataName);
end
disp('Setting up Shaded Plots')

%% Here we go

% Define Time frame of interest

TOF= [-0.1 0.35];

Time= dummy.time;
[TimeInd]= find(Time >= TOF(1) & Time <= TOF(2));
Time= Time(TimeInd);

% Define data of interest

Labels{1}= {'STDs', 'STDw'};
Labels{2}= {'DEVs', 'DEVw'};
Labels{3}= {'STD', 'DEV'};

Contrasts = {'STDs-w', 'DEVs-w', 'STD-DEV'};


%% Single-subject plots

fprintf('\nPreparing Single-subj plots')

for ll = 1:length(Labels)
    
    figure('Color', 'White'); hold on
    
    for ee= 1:length(Labels{ll})
        
        data= getfield(dummy, Labels{ll}{ee});
        
        if ee == 1, Colors= '-b';    else, Colors= '-r';    end
        
        n = 1;
        for ss = 1:2:NSubj
            
            subplot(4,5,n); hold on
            data1.time = dummy.time;
            data1.nanmean = nanmean(data(ss:ss+1,:));
            
            %% Let's plot
            
            Minvals= [];
            Maxvals= [];
            
            % Single plots
            
            plot(data1.time, data1.nanmean, Colors);
            title([Contrasts{ll} sprintf('SubjN%d', n)])
            n = n+1;
            %     legend(Labels)
            
            if ee == 2
                % Add lines
                plot1=gcf; plot1= plot1.CurrentAxes;
                plot1.FontSize = 7;
                plot1.FontWeight = 'bold';
                
                Minvals= min(plot1.YLim);
                Maxvals= max(plot1.YLim);
                
                ylimit = [min(Minvals) max(Maxvals)];
                %[min(Minvals)+0.01*min(Minvals) max(Maxvals)+0.01*max(Maxvals)];
                ylim(ylimit);
                
                hline = refline([0,0]);
                hline.Color = 'k';
                line([0,0], ylimit, 'color','black','LineStyle','--');
                line([0.3,0.3], ylimit, 'color','black','LineStyle','--');
                
                xlabel('\fontsize{9} Time (s)', 'fontweight', 'bold')
                ylabel('\fontsize{9} Amplitude uV', 'fontweight', 'bold')
                
            end
            
        end
    end
end


%% Plotting GAs

fprintf('\nPreapring GA plots')

for ll = 1:length(Labels)
    
    figure('Color', 'White'); hold on
    
    for ee= 1:length(Labels{ll})
        
        data= getfield(dummy, Labels{ll}{ee});
        
        if ee == 1, Colors= '-b';    else, Colors= '-r';    end
        
        data1= [];
        data1.time= Time;
        data1.nanmean= nanmean(data);
        data1.SE= (std(data))/sqrt(NSubj);
        
        %% Let's plot
        
        Minvals= [];
        Maxvals= [];
        
        % Single plots
        
        h{ee} = shadedErrorBar(data1.time, data1.nanmean, data1.SE, 'lineprops', Colors);
        title(Contrasts{ll})
        %     legend(Labels)
        
        if ee == 2
            
            leg = legend([h{1}.patch, h{2}.patch],Labels{ll});
            leg.AutoUpdate = 'off';
            leg.Box = 'off';
            
            % Add lines
            plot1=gcf; plot1= plot1.CurrentAxes;
            plot1.FontSize = 12;
            plot1.FontWeight = 'bold';
            
            Minvals= min(plot1.YLim);
            Maxvals= max(plot1.YLim);
            
            ylimit = [min(Minvals) max(Maxvals)];
            %[min(Minvals)+0.01*min(Minvals) max(Maxvals)+0.01*max(Maxvals)];
            ylim(ylimit);
            
            hline = refline([0,0]);
            hline.Color = 'k';
            line([0,0], ylimit, 'color','black','LineStyle','--');
            line([0.3,0.3], ylimit, 'color','black','LineStyle','--');
            
            xlabel('\fontsize{12} Time (s)', 'fontweight', 'bold')
            ylabel('\fontsize{12} Amplitude uV', 'fontweight', 'bold')
            
        else
            Data1 = data;
        end
        
    end
    
        % Statistical testing
        
        alpha = .05; %alpha/nperm
        [H,P,CI,STATS] = ttest(Data1, data);

        % Adjust pvals
        [h1, p_crit, adj_ci_cvrg, adj_p]= fdr_bh_pCorr(P,.05,'pdep');
      
        % Plot stats
        Stars = nan(size(adj_p));
        Stars(find(adj_p < p_crit)) = 1.8;
        plot(data1.time, Stars, '*k');
        leg.String{3} = 'pAdj <.05';
        
        clear p observeddifference effectsize p2Adjust Data1 data
        
end

fig2save = 'C:\Users\anton\OneDrive\Work\UM\BAND\Writing\Subjective Rhythmization\Humans - BeatIt\Submission\Figures\ERPs';
saveFigures_AC(fig2save)

