%% 统计特征 | statistical features
%% 最小值 | Minimum
% X: single channel EEG signal (either a row vector or a column vector)


function X_min = feat_Minimum(X,~)
X_min = min(X);
end

