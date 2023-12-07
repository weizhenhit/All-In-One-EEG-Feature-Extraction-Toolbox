%% 统计特征 | statistical features
%% 二阶差分平均值 | Mean of Second Difference
% 描述了时域信号强度的变化
% X: single channel EEG signal (either a row vector or a column vector)

%% Reference
%       [1] Jenke, R., Peer, A., and Buss, M. (2014). Feature Extraction and Selection for Emotion Recognition from EEG. IEEE Transactions on Affective Computing 5, 327–339. doi: 10.1109/TAFFC.2014.2339834.

function SD = feat_SecondDifference(X,~)
T = length(X); 
Y = 0;
for t = 1 : T - 2
  Y = Y + abs(X(t+2) - X(t));
end
SD = (1 / (T - 2)) * Y;
end