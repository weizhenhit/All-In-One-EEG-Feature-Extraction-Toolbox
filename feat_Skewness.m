%% （高阶）统计特征 | statistical features
%% 偏度 | Skewness
% 表征信号在质心附近的不对称性
% The use of these features is motivated by the fact that distribution of EEG signal is often characterized by its level of dispersion, asymmetry, and concentration around the mean.
% X: single channel EEG signal (either a row vector or a column vector)

function SKEW = feat_Skewness(X,~)
SKEW = skewness(X);
end

