
%% Processing lines

% Individual neurophysiological signatures of spontaneous rhythm processing
% A. Criscuolo, M. Schwartze, M.J. Henry, C. Obermeier & S.A. Kotz
% NeuroImage, 2023
% https://doi.org/10.1016/j.neuroimage.2023.120090


% 1. Time-locked ERP analyses 

Preprocess_AC(1, 1, 2)
Process_AC(1, 1, 0, 2)

% Plot event-locked ERP data
ShadedPlots_ERP_CompleteData(1, 0) % sess1, no BSL corr


% 2. Time-locked TFR analyses

Preprocess_AC(1, 2, 2) 
Process_AC(1, 2, 0, 2)

% Plot event-locked TFR data
ShadedPlots_TFR_CompleteData(1, 0)
PowerPlots_StemPlots_AC(1, 0)


% 3. Beat (accent) modelling based on TFR data
modelBeat_AC(1, 0) % includes plots and group-level stats

% DEV plots according to Beat modelling
Plots_ERP_BinaryBeat(1, 0) 


% 4. Analyses on continuous data: 
% - FFT; 
% - phase analyses.

Data4FFT_AC(1) % segment over the entire auditory sequence and use fourier transform
EntrainME_AC(1) % extract signal at stimulus periodicity and perform phase-locking analyses


% More figures in
% figures4slide

