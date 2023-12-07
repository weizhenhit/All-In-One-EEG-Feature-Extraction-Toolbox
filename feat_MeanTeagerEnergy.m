%% 非线性特征 | nonlinear features 
%% Mean teager energy (MTE)
% X: single channel EEG signal (either a row vector or a column vector)
% See also:
%       feat_ShannonEntropy, feat_LogEnergyEntropy, feat_TsallisEntropy, feat_RenyiEntropy, feat_MeanEnergy

%% Reference
%       [1] Şen, B., Peker, M., Çavuşoğlu, A., and Çelebi, F. V. (2014). A Comparative Study on Classification of Sleep Stage Based on EEG Signals Using Feature Selection and Classification Algorithms. J Med Syst 38, 18. doi: 10.1007/s10916-014-0018-0.
%       [2] Gardner, A. B., Krieger, A. E., Vachtsevanos, G., and Litt, B., Oneclass novelty detection for seizure analysis from intracranial EEG. J Mach Learn Res 7:1025–1044, 2006

function MTE = feat_MeanTeagerEnergy(X,~)
N = length(X); 
Y = 0;
for m = 3:N
  Y = Y + ((X(m-1) ^ 2) - X(m) * X(m-2));
end
MTE = (1 / N) * Y;
end

