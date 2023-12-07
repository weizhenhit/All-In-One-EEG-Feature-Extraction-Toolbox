%% 非线性（动力学）特征：基于时域的熵 | nonlinear features: time-domain-based entropy 
%% Renyi Entropy
% X: single channel EEG signal (either a row vector or a column vector)
% See also:
%       feat_ShannonEntropy, feat_LogEnergyEntropy, feat_TsallisEntropy, feat_MeanEnergy, feat_MeanTeagerEnergy

%% Reference
%       [1] Şen, B., Peker, M., Çavuşoğlu, A., and Çelebi, F. V. (2014). A Comparative Study on Classification of Sleep Stage Based on EEG Signals Using Feature Selection and Classification Algorithms. J Med Syst 38, 18. doi: 10.1007/s10916-014-0018-0.
%       Rényi entropy, introduced by Rényi [50], is a special case of spectral entropy based on the concept of generalized entropy of a probability distribution. If p is a probability distribution on a finite set, its RE of order q is defined to be RE is calculated as q=2 in this study.
%       [2] Hariharan, M., Sindhu, R., Vijean, V., Yazid, H., Nadarajaw, T., Yaacob, S., et al. (2018). Improved binary dragonfly optimization algorithm and wavelet packet based non-linear features for infant cry classification. Comput Meth Prog Bio 155, 39–51. doi: 10.1016/j.cmpb.2017.11.021.

function ReEn = feat_RenyiEntropy(X,opts) 

    if isfield(opts,'alpha')	% 指定
        alpha = opts.alpha; 
    else	% 未指定 default
        alpha = 2;
    end
    % Convert probability using energy
    P = (X .^ 2) ./ sum(X .^ 2);     % 计算每个样本的概率
    % Entropy 
    En = P .^ alpha; 
    ReEn = (1 / (1 - alpha)) * log2(sum(En)); 
end

