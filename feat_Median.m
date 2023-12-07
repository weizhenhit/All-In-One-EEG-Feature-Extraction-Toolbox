%% 统计特征 | statistical features
%% 中位数 | Median
% X: single channel EEG signal (either a row vector or a column vector)

function X_med = feat_Median(X,~)
X_med = median(X);
end

