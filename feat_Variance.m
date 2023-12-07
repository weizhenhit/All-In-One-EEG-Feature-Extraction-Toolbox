%% 统计特征 | statistical features
%% 方差 | Variance
% X: single channel EEG signal (either a row vector or a column vector)

function VAR = feat_Variance(X,~)
VAR = var(X);
end
