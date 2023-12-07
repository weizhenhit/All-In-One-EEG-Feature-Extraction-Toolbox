%% 统计特征 | statistical features
%% 最大值 | Maximum
% X: single channel EEG signal (either a row vector or a column vector)


function X_max = feat_Maximum(X,~)
X_max = max(X);
end

