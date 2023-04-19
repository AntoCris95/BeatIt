

% This function runs ICA by using the FieldTrip toolbox %

% The function is called during Preprocess_AC.
% There are an automatic and manual version. Specify the preference at the
% beginning of the function.
% In the manual version, you can inspect the correlation between the
% time-course of each component and the EEG data and manually inputs the
% components that you want to remove. Moreover, you can inspect the
% topography of correlations for a further check.
% In the automatic version, instead, the algorithm will automatically
% select components who have a correlation above a given specified
% threshold and would remove them. However, in case there are more that a
% give N of components as indicated, ICA is run over again.
% Data and components to remove are saved at the end.


% Written in July 2020
% Created by Antonio Criscuolo



function [data] = goICA_AC(data, outdir, Name2save, Auto)

if ~exist('Auto', 'var')
    Auto = 0; % automatic ICA = 1; manual = 0;
end

Fullname2save = fullfile(outdir, Name2save);
Comp2save = fullfile(outdir, ['Components_' Name2save]);

if ~exist(Comp2save)
    
    % Preprocess channels of interest
    
    cfg = [];
    cfg.channel = {'EOG'}; % select EOG ch
    cfg.continuous = 'yes';
    dEOG = ft_preprocessing(cfg, data);
    
    cfg = [];
    cfg.channel = {'EEG'}; % select EEG ch
    cfg.continuous = 'yes';
    data = ft_preprocessing(cfg, data);
    
    % RUN ICA on distributed computing
    
    cfg            = [];
    cfg.channel    = {'EEG'};
    cfg.method     = 'fastica'; % 'runica'
    cfg.demean     = 'no';
    cfg.numcomponent = 20; % or 'all'
    comp           = ft_componentanalysis(cfg, data);
    
    %Save components
    
    save(Comp2save,'dEOG', 'comp');
    
else
    load(Comp2save)
end

% Time-locked components

% Here we run time-locked analysis on the above specified channels.
% Then we calculate the correlation between the time-course
% of the signals and plot it in a histogram.
% Based on correlations, we manually select a number of components
% exceeding >.4 correlation and inspect topographies.
% Based on the topographic location of these selected components,
% we select the ones to reject.

%EOG
cfg = [];
cfg.keeptrials = 'yes';
tlEOG = ft_timelockanalysis(cfg, dEOG);
tlComp = ft_timelockanalysis(cfg, comp);

x = tlEOG.trial(:,1,:); % blinks
for c = 1:size(tlComp.trial,2)
    y = tlComp.trial(:,c,:); % components
    rBlink(c) = corr(y(:), x(:));
end

if Auto == 0
    
    figure;
    bar(1:c, abs(rBlink), 'b'); title('BLINKS');
    xlabel('comp #');
    
    % Note down the components to inspect
    % NB: make a vector [].
    
    comp2rem = input('bad components are: ')
    
else
    
    comp2rem = find(abs(rBlink) > 0.35);
    
    if length(comp2rem) >= 3 || length(comp2rem) < 1
        
        % Re-run ICA
        cfg            = [];
        cfg.method     = 'fastica';
        cfg.demean     = 'no';
        cfg.numcomponent = 20; % or 'all'
        comp           = ft_componentanalysis(cfg, data);
        
        %EOG
        cfg = [];
        cfg.keeptrials = 'yes';
        tlEOG = ft_timelockanalysis(cfg, dEOG);
        tlComp = ft_timelockanalysis(cfg, comp);
        
        x = tlEOG.trial(:,1,:); % blinks
        for c = 1:size(tlComp.trial,2)
            y = tlComp.trial(:,c,:); % components
            rBlink(c) = corr(y(:), x(:));
        end
        
        comp2rem = find(abs(rBlink) > 0.35);
        
    end
end

% Look at components' topographies

if Auto == 0
    if comp2rem > 0
        cfg           = [];
        cfg.component = [comp2rem];
        cfg.layout    = 'EEG1020.lay';
        cfg.comment   = 'no';
        figure; ft_topoplotIC(cfg, comp)
        
        comp2rem = input('bad components are: ')
    end
end

% Save components

save(Comp2save,'dEOG', 'comp', 'comp2rem');

% Reject components

if comp2rem > 0
    cfg           = [];
    cfg.component = comp2rem;
    data = ft_rejectcomponent(cfg, comp,data);
end

% Remove eog channels
selchan = ft_channelselection({'all', '-EOG'}, data.label);
cfg = [];
cfg.channel = selchan;
data = ft_selectdata(cfg, data);

% Save ICA data

if ~isreal(data.trial{1})
    disp('ICA could not be performed. Interrupting')
    return
else
    save(fullfile(outdir, Name2save), 'comp2rem','data');
end
