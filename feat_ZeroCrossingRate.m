%% 统计特征 | statistical features
%% 过零率 | Zero-crossing rate
% X: single channel EEG signal (either a row vector or a column vector)

%% Reference
%       [1] Şen, B., Peker, M., Çavuşoğlu, A., and Çelebi, F. V. (2014). A Comparative Study on Classification of Sleep Stage Based on EEG Signals Using Feature Selection and Classification Algorithms. J Med Syst 38, 18. doi: 10.1007/s10916-014-0018-0.


function zero_crossing = feat_ZeroCrossingRate(X,~)
    T = length(X);
    zero_crossing = 0;
    for i = 1:T-1
        if(X(i) * X(i+1) < 0)   % 只要向量两个样本点的符号相反，zero-crossing++
            zero_crossing = zero_crossing + 1;
        end
    end
    zero_crossing = zero_crossing / (T-1);
end

% 更简便的实现方式
% function zero_crossing = feat_ZeroCrossingRate2(X,~)
%     T = length(X);
%     zero_crossing = sum(abs(diff(X > 0))) / (T-1);
% end