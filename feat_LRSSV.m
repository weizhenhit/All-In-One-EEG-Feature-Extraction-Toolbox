%% 统计特征 | statistical features
%% Log Root Sum Of Sequential Variation
% The log root sum of sequential variations (LRSSV) is proposed in this paper to measure the sequential variations between the samples of the signal.
% X: single channel EEG signal (either a row vector or a column vector)

%% Reference
%       [1] Memar, P., and Faradji, F. (2018). A Novel Multi-Class EEG-Based Sleep Stage Classification System. IEEE Transactions on Neural Systems and Rehabilitation Engineering 26, 84–95. doi: 10.1109/TNSRE.2017.2776149.

function LRSSV = feat_LRSSV(X,~)
N = length(X); 
Y = zeros(1, N-1);
for i = 2:N
  Y(i-1) = (X(i) - X(i-1)) ^ 2;
end
LRSSV = log10(sqrt(sum(Y)));
end
