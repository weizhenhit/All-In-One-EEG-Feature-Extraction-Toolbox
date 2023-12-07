%% 微分熵 | differential entropy (DE)
% 微分熵实际上是香农熵在连续信号上的推广，微分熵并不像上文单指一种频带分段方法，而是对特征附加的非线性变换，尤其在频域能量特征中效果显著
% 简化计算：直接对能量特征求log
% X: single channel EEG signal (either a row vector or a column vector)

%% Reference
%       [1] https://blog.csdn.net/zhoudapeng01/article/details/108596687 (theory)
%       [2] https://blog.csdn.net/YuqingF/article/details/125873032 (code)


function DE = feat_DifferentialEntropy(X, ~)
    variance = var(X);   % 默认情况下，方差按观测值数量N -1 实现归一化
    DE = log(2 * pi * exp(1) * variance) / 2;    % 简化版的微分熵求解公式
end
