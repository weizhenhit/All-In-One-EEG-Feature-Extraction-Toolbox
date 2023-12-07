%% 非线性（动力学）特征：基于时域的熵 | nonlinear features: time-domain-based entropy 
%% Log Energy Entropy
% 表征时间序列数据的复杂性（complexity）In some studies, the log energy entropy is used as a feature to distinguish focal(病灶的) and non-focal EEG signals.
% X: single channel EEG signal (either a row vector or a column vector)
% See also:
%       feat_ShannonEntropy, feat_RenyiEntropy, feat_TsallisEntropy,
%       feat_MeanEnergy, feat_MeanTeagerEnergy

%% Reference
%       [1] Gupta, V., Priya, T., Yadav, A. K., Pachori, R. B., and Rajendra Acharya, U. (2017). Automated detection of focal EEG signals using features extracted from flexible analytic wavelet transform. Pattern Recognition Letters 94, 180–188. doi: 10.1016/j.patrec.2017.03.017.

function LogEn = feat_LogEnergyEntropy(X,~)
    LogEn = sum(log(X .^ 2)); 
end

