function cfg = config_analysis()
% CONFIG_ANALYSIS returns all configurable parameters for this project.

    % 
    cfg.projectRoot = fileparts(mfilename('fullpath'));
    cfg.dataRoot    = fullfile(cfg.projectRoot, 'data', 'raw');
    cfg.derivRoot   = fullfile(cfg.projectRoot, 'derivatives');
    cfg.resultsRoot = fullfile(cfg.projectRoot, 'results');
    cfg.figureDir   = fullfile(cfg.resultsRoot, 'figures');

    % 
    % readtable(fullfile(cfg.projectRoot,'patient_groups.csv'))
    cfg.groups.control = { 'sub-HC-001', 'sub-HC-002' ...};   
    cfg.groups.good    = { 'sub-MCS-001', 'sub-MCS-002' ...}; % 
    cfg.groups.bad     = { 'sub-UWS-001', 'sub-UWS-002' ...}; % 

    % 
    cfg.nChannels           = 19;
    cfg.analysisDurationSec = 5 * 60;   
    cfg.freqBands = [1 4; 4 8; 8 13; 13 30; 30 50; 50 70];
    cfg.winLenSec = 4;
    cfg.stepSec   = 1;

    % 
    cfg.nClusters = 4;  

    % 
    cfg.randomSeed = 42;

    % 
    cfg.groupLabels  = {'HC', 'MCS', 'UWS'};
    cfg.featureNames = {'GlobalVariability','Metastability', ...
                        'ComplexityMean','ComplexityStd'};
end