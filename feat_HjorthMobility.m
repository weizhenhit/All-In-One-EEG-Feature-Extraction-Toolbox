%% Hjorth parameters
% Hjorth Activity 活动性：衡量信号波幅的偏离程度，是时间函数的方差var
% Hjorth Mobility 移动性：衡量坡度的变化，为信号一阶导数方差的平方根 除以 信号方差的平方根
% Hjorth Complexity 复杂性：衡量一个振幅上有多少个标准的坡
% X: single channel EEG signal (either a row vector or a column vector)
% See also:
%       feat_HjorthActivity, feat_HjorthComplexity

%% Reference
%       [1] Şen, B., Peker, M., Çavuşoğlu, A., and Çelebi, F. V. (2014). A Comparative Study on Classification of Sleep Stage Based on EEG Signals Using Feature Selection and Classification Algorithms. J Med Syst 38, 18. doi: 10.1007/s10916-014-0018-0.
%       [2] Jenke, R., Peer, A., and Buss, M. (2014). Feature Extraction and Selection for Emotion Recognition from EEG. IEEE Transactions on Affective Computing 5, 327–339. doi: 10.1109/TAFFC.2014.2339834.


function HM = feat_HjorthMobility(X,~)
    % First derivative
    x0  = X(:);
    x1  = diff([0; x0]); 
    % Standard deviation 
    sd0 = std(x0); 
    sd1 = std(x1);
    HM  = sd1 / sd0;
end

