%% 统计特征 | statistical features
%% 标准差 | Standard Deviation
% X: single channel EEG signal (either a row vector or a column vector)

function SD = feat_StandardDeviation(X,~)
N  = length(X); 
mu = mean(X); 
SD = sqrt((1 / (N - 1)) * sum((X - mu) .^ 2));
% SD = std(X);
end

