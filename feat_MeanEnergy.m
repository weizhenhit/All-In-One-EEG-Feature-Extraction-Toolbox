%% 非线性特征 | nonlinear features 
%% 平均瞬时能量 Mean Energy
% X: single channel EEG/imfx/hilbert(imfx) signal (either a row vector or a column vector)
% 对于signal或者imf这样的时域信号，abs(signal).^2得到的是每个点处的瞬时能量（即瞬时振幅的平方，是一个长度为N的向量）,再取mean后得到的是平均瞬时能量（是一个值）
% 对于hilbert(imf1)这样得到的是imf的解析信号表示（是个复数），abs(hilbert(imf1)).^2得到的是每个点处的瞬时能量（同样是瞬时振幅的平方，是一个长度为N的向量），再取mean后得到的是平均瞬时能量（是一个值）
% See also:
%       feat_ShannonEntropy, feat_LogEnergyEntropy, feat_TsallisEntropy, feat_RenyiEntropy, feat_MeanTeagerEnergy

%% Reference
%       [1] Şen, B., Peker, M., Çavuşoğlu, A., and Çelebi, F. V. (2014). A Comparative Study on Classification of Sleep Stage Based on EEG Signals Using Feature Selection and Classification Algorithms. J Med Syst 38, 18. doi: 10.1007/s10916-014-0018-0.
%       Mean energy (ME) Energy values in signals increase along with increases in activity. When different activities in different sleep stages are considered, it is believed that the mean energy is a good indicator. It is calculated as the average of the squares of all samples in the signal.

function ME = feat_MeanEnergy(X,~)
    ME = mean(abs(X) .^ 2);
end