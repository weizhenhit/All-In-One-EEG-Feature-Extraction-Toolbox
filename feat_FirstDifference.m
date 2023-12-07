%% 统计特征 | statistical features
%% 一阶差分平均值 | Mean of First Difference/Mean Curve Length（MCL）
% 描述了时域信号强度的变化，是一种量化信号复杂性的简单度量，在分析和经验上与Katz分形维数相关，上述两者是同一个东西，可用于区分不同被试
% X: single channel EEG signal (either a row vector or a column vector)

%% Reference
%       [1] Yahyaei, R., and Esat Özkurt, T. (2022). Mean curve length: An efficient feature for brainwave biometrics. Biomedical Signal Processing and Control 76, 103664. doi: 10.1016/j.bspc.2022.103664.


function FD = feat_FirstDifference(X,~)
T = length(X);
Y = 0;
for t = 1 : T - 1
  Y = Y + abs(X(t+1) - X(t));
end
FD = (1 / (T - 1)) * Y;
end