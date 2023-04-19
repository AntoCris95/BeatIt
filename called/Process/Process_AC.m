
%% Processing function for both ERP and TFR

% Call the function after data has been preprocessed.
% This processing function performs:

% 1. Bsl/Mean/No correction;
% 2. Average epochs;
% 3. Grand-average;

% In case of TFR data, however, data first undergoes frequency-transform and
% then is fed into the pipeline. 

% Consult the analysis pipeline for info.

% Created in March 2020
% Edited in March 2021
% Written by Antonio Criscuolo


function Process_AC(sess, ERPTFR, datatype, BSLcorr)

if nargin < 1
    sess = input('Session number?');
    ERPTFR = input('ERP (1) or TFR (2) analysis?');
    if ERPTFR == 2
    datatype = input('Are these simple (0) or complex data (1)?');
    else, datatype = 0;
    end
    BSLcorr = input('Apply BSL correction (1), mean correction (2) or nocorr (0)?');
end

% Here we go

disp('Setting up Preprocessing')

if ERPTFR == 2
    
    disp('Preparing TF-analyses')
    go_FastTFR_AC(sess)
%     TFR_Fieldtrip(sess, datatype)
    
    disp('Performing Bsl/Mean correction and epochs averaging')
    MeanCorr_AvgEp_TFR_AC(sess, datatype, BSLcorr)

else

    disp('Performing Bsl/Mean correction and epochs averaging')
    MeanCorr_AvgEp_ERP_AC(sess, BSLcorr)

end

disp('Moving to GA')
GA(sess, ERPTFR, datatype, BSLcorr)
GGA(sess, ERPTFR, datatype, BSLcorr)


% Emit sound when the pipeline is completed

% sound(cos(1:5000))
if ERPTFR == 2
    disp('Might need hp filter 2Hz before stats')
end

