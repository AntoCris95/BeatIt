
function data = get_TFR_AC(fullfilename, outdir, datatype, FQS)

% Set up filename and paths
datadir = strsplit(fullfilename, 'ep'); datadir = datadir{1};
filename = erase(fullfilename, datadir);

load(fullfile(fullfilename));

if datatype % complex exponentials
    saveName= fullfile(outdir, ['cwtcomplex ' filename]);
else % real numbers only
    saveName= fullfile(outdir, ['cwt ' filename]);
end

% Set up params

if ~exist('FQS', 'var')
    FQStep = 1;
    FQS = 1 : FQStep : 40;
end
if length(FQS) == 2
    FQStep = 1;
    FQS = FQS(1) : FQStep : FQS(2);
end

Method = 'wavelet'; % options: mtmconvol, tfr, wavelet

windowsize = data.fsample/10e3; % e.g., 25ms
TOI = data.time{1}(1) + (windowsize) : windowsize : data.time{1}(end) - (windowsize);

if ~exist(saveName)

    % Procede to time-frequency analysis

    cfg = [];
    cfg.method = Method;

    if datatype
        cfg.output = 'fourier'; % complex spectra
    else
        cfg.output = 'pow';
    end

    % Standard params
    cfg.channel         = {'EEG'};
    cfg.keeptrials      = 'yes';
    cfg.pad             = 'nextpow2';
    cfg.polyremoval     = 0; % the 0 order polyn (mean subtraction; the 1st order (detrending + demean)
    cfg.toi             = TOI;
    cfg.foi             = FQS;

    % Call distributed computing
    cfg.inputfile       = fullfilename;
    cfg.outputfile      = saveName; % save

    % Method-specific params

    if strcmp(cfg.method, 'mtmconvol')

        smoothing       = 2; % only effective
        cfg.taper       = 'hanning';
        cfg.tapsmofrq   = smoothing;
        cfg.t_ftimwin   = ones(length(cfg.foi),1); % 1x numfoi - length of the time-window, in seconds

        data = ft_freqanalysis(cfg);

    else % continuous wavelet transform

        data = call2wavelet_AC(cfg);

    end

else
    load(saveName)
    data = freq; clear freq
end

