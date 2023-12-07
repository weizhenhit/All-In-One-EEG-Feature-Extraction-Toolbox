%% 非线性（动力学）特征：基于时域的熵 | nonlinear features: time-domain-based entropy 
%% 香农熵 | Shannon Entropy 
% 表征时间序列数据的复杂性（complexity）和规律性（regularity），值越大，时间序列越复杂、越不规律
% X: single channel EEG signal (either a row vector or a column vector)
% See also:
%       feat_ApproximateEntropy, feat_SampleEntropy,
%       feat_PermutationEntropy, feat_FuzzyEntropy, feat_LogEnergyEntropy,
%       feat_RenyiEntropy, feat_TsallisEntropy, feat_MeanEnergy, feat_MeanTeagerEnergy

function ShEn = feat_ShannonEntropy(X,~) 
    % Convert probability using energy
    P    = (X .^ 2) ./ sum(X .^ 2);     % 计算每个样本的概率
    % Entropy 
    En   = P .* log2(P);
    ShEn = -sum(En); 
end

