function run_dFC_preprocessing(cfg)
% RUN_DFC_PREPROCESSING calculates the dynamic wPLI for each subject and saves the results.
% Usage：cfg = config_analysis(); run_dFC_preprocessing(cfg);

    if nargin < 1
        cfg = config_analysis();
    end
    rng(cfg.randomSeed);

    if ~exist(cfg.derivRoot, 'dir'), mkdir(cfg.derivRoot); end

    allPatients = [cfg.groups.control, cfg.groups.good, cfg.groups.bad];

    for pa = 1:numel(allPatients)
        patient = allPatients{pa};

        inputFile = fullfile(cfg.dataRoot, patient, 'process', 'EEG.mat');
        if ~isfile(inputFile)
            warning('File not found. Skipping. ：%s', inputFile);
            continue;
        end

        S = load(inputFile, 'EEG');
        if ~isfield(S, 'EEG')
            warning('No EEG variables were found, so skipping. ：%s', inputFile);
            continue;
        end
        EEG = S.EEG;

        fs = EEG.srate;
        nSamples = min(size(EEG.data, 2), cfg.analysisDurationSec * fs);
        data_normal = EEG.data(1:cfg.nChannels, 1:nSamples);

        nBands = size(cfg.freqBands, 1);
        data_freq = cell(1, nBands);
        for fb = 1:nBands
            data_freq{fb} = eegfilt(data_normal, fs, ...
                cfg.freqBands(fb,1), cfg.freqBands(fb,2));
        end

        winSamples  = round(cfg.winLenSec * fs);
        stepSamples = round(cfg.stepSec * fs);
        start_idx   = 1 : stepSamples : (nSamples - winSamples + 1);
        nWindows    = numel(start_idx);

        wpli_dynamic = cell(1, nBands);
        for fre = 1:nBands
            wpli_dynamic1 = zeros(cfg.nChannels, cfg.nChannels, nWindows);
            for nw = 1:nWindows
                idx = start_idx(nw):(start_idx(nw)+winSamples-1);
                data1 = data_freq{fre}(:, idx);
                data2 = zscore(data1')';
                data3 = reshape(data2, [size(data2,1), 1, size(data2,2)]);

                [wpli1, ~] = conn_wPLI_hao(data3);
                wpli1(isnan(wpli1)) = 0;
                wpli2 = squeeze(wpli1);
                wpli3 = wpli2 + wpli2';
                wpli_dynamic1(:,:,nw) = wpli3;
            end
            wpli_dynamic{fre} = wpli_dynamic1;
        end

        outputDir  = fullfile(cfg.derivRoot, patient);
        if ~exist(outputDir, 'dir'), mkdir(outputDir); end
        outputFile = fullfile(outputDir, 'wpli_dynamic.mat');
        save(outputFile, 'wpli_dynamic', '-v7.3');
    end
end