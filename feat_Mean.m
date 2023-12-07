%% 统计特征 | statistical features
%% 平均值 | Mean
% X: single channel EEG signal (either a row vector or a column vector)


function Mean = feat_Mean(X,~)
Mean = mean(X);
end

