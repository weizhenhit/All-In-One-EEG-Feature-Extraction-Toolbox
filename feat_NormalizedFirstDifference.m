%% 统计特征 | statistical features
%% 归一化一阶差分 | Normalized First Difference
% X: single channel EEG signal (either a row vector or a column vector)

%% Reference
%       [1] Jenke, R., Peer, A., and Buss, M. (2014). Feature Extraction and Selection for Emotion Recognition from EEG. IEEE Transactions on Affective Computing 5, 327–339. doi: 10.1109/TAFFC.2014.2339834.

function NFD=feat_NormalizedFirstDifference(X,~)
T = length(X); 
Y = 0; 
for t = 1 : T - 1
  Y = Y + abs(X(t+1) - X(t));
end
FD  = (1 / (T - 1)) * Y; 
NFD = FD / std(X);
end