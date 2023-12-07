clear all

load('X.mat');      % single channel EEG signal (either a row vector or a column vector)
% X = X';

%% You can call and test each feature extraction function separately.
% power spectral density (PSD), welch method
opts.fs = 128;      % sampling frequency
opts.f_low = 1; opts.f_high = 45;   % frequency range
PSD = feat_PSDWelch(X, opts);   disp(PSD);
% empirical mode decomposition (EMD)
opts.imf_level = 3;     % extracted maximum number of intrinsic mode functions
imf_features = feat_EmpiricalModeDecomposition(X, opts);    disp(imf_features);

%% You can also use the 'feat_main' function to extract multiple features at once.
% The idea of feat_main.m comes from https://github.com/JingweiToo/EEG-Feature-Extraction-Toolbox
opts.fs = 128;      % sampling frequency
opts.f_low = 1; opts.f_high = 45;   % feat_PSDPer, feat_PSDWelch, feat_PSDArburg
opts.imf_level = 3;         % fet_EmpiricalModeDecomposition
opts.f_bin = [1 45];        % fet_HilbertHuangSpectralEntropy
opts.wavefunc = 'db4';      % feat_WaveletEntropy

feature_abbr = {'mean','sd','var','md','max','min','1d','n1d','2d','n2d','kurt','skew','c0r','lrssv','ha','hm','hc','sh','se','ae','pe','fe','re','te','le','me','mte','fd','plzc','psd_welch','de','emd','hhse','we'}; % 34
feature_set = [];   % 104
for feat = 1:length(feature_abbr)
    disp(feature_abbr{feat});
    tic
    feature_set = [feature_set feat_main(feature_abbr{feat}, X, opts)];
    toc
end

%% Common Spatial Patterns (CSP): feat_MulticlassCSP1v1_seg
load('EEG.mat')   % data from SEED datasset, EEG_train and EEG_test, times=1000, channels=62, nbClasses=3
%   EEG_train: 训练集各类别的所有trials：
%       EEG_train.x: EEG，times*channels*trials
%       EEG_train.y: 训练集所有样本的类别标签，trials*1
%   EEG_test: 测试集各类别的所有trials：
%       EEG_test.x: EEG，times*channels*trials
%       EEG_test.y: 测试集所有样本的类别标签，trials*1
opts.fs = 200;
opts.m = 4;     % number of CSP filters, dimensions of CSP features are 2m, m<channels
opts.seg_len = 1;   % length of each segments，sec
[feature_train, feature_test, CSPMatrix] = feat_MulticlassCSP1v1_seg(EEG_train, EEG_test, opts);
% feature_train: (trials, m*nbClasses*(nbClasses-1), nsegments)