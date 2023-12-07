%% （高阶）统计特征 | statistical features
%% 峰度 | kurtosis
% 表征信号在质心附近的平坦度
% The use of these features is motivated by the fact that distribution of EEG signal is often characterized by its level of dispersion, asymmetry, and concentration around the mean.
% X: single channel EEG signal (either a row vector or a column vector)


function KURT = feat_Kurtosis(X,~)
% Kurtosis 
KURT = kurtosis(X);
end
