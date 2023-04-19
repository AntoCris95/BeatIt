
function data = call2wavelet_AC(cfg)

FQS         = cfg.foi;
FQStep      = FQS(2)-FQS(1);

if FQS(end) > 5 % use different N of cycles for low and high freqs

    cfg.foi     = [FQS(1):FQStep:5];
    cfg.width   = 3; % low-frequency with 3 cycles

    data1        = ft_freqanalysis(cfg);

    cfg.foi     = [6:FQStep:FQS(2)];
    cfg.width   = 7; % high-frequency with 7 cycles

    data2       = ft_freqanalysis(cfg);

    % Concatenate the 2 datasets

    if any(ismember(fieldnames(data1), 'powspctrm'))
        datatype = 0;
    else
        datatype = 1;
    end

    if datatype

        Nfqs= size(data2.fourierspctrm, 3);
        data1.fourierspctrm = data1.fourierspctrm;
        data1.fourierspctrm(:,:,end+1 : end+Nfqs, :) = data2.fourierspctrm;

    else

        Nfqs= size(data2.powspctrm, 3);
        data1.powspctrm = data1.powspctrm;
        data1.powspctrm(:,:,end+1 : end+Nfqs, :) = data2.powspctrm;

    end

    data = data1; % to send back

else
    cfg.width   = 3; % low-frequency with 3 cycles
    data = ft_freqanalysis(cfg);
end

data.freq = FQS;
