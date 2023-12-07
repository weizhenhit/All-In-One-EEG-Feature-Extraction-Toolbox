%% 非线性（动力学）特征：基于时域的熵 | nonlinear features: time-domain-based entropy 
%% Tsallis Entropy
% X: single channel EEG signal (either a row vector or a column vector)
% See also:
%       feat_ShannonEntropy, feat_LogEnergyEntropy, feat_RenyiEntropy,
%       feat_MeanEnergy, feat_MeanTeagerEnergy

%% Reference
%       [1] Hariharan, M., Sindhu, R., Vijean, V., Yazid, H., Nadarajaw, T., Yaacob, S., et al. (2018). Improved binary dragonfly optimization algorithm and wavelet packet based non-linear features for infant cry classification. Comput Meth Prog Bio 155, 39–51. doi: 10.1016/j.cmpb.2017.11.021.

function TsEn = feat_TsallisEntropy(X,opts)
    if isfield(opts,'alpha')	% 指定
        alpha = opts.alpha; 
    else	% 未指定 default
        alpha = 2;
    end
    % Convert probability using energy 
    P = (X .^ 2) ./ sum(X .^ 2);     % 计算每个样本的概率
    % Entropy 
    En = P .^ alpha;
    TsEn = (1 / (alpha - 1)) * (1 - sum(En(:))); 
end

