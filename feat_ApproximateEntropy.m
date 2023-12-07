%% 非线性（动力学）特征：基于时域的熵 | nonlinear features: time-domain-based entropy 
%% 近似熵 | Approximate Entropy （similar to Sample Entropy）
% 表征时间序列数据的复杂性（complexity）和规律性（regularity），值越大，时间序列越复杂、越不规律
% X: single channel EEG signal (either a row vector or a column vector)
% opts.r：相似容限阈值，通常r在0.1std(X)和0.25std(X)之间时，ApEn有比较合理的统计特性，默认0.2std(X)
% See also:
%       feat_ShannonEntropy, feat_SampleEntropy, feat_PermutationEntropy, feat_FuzzyEntropy

%% Reference
%       [1] Kijoon Lee (2023). Fast Approximate Entropy (https://www.mathworks.com/matlabcentral/fileexchange/32427-fast-approximate-entropy), MATLAB Central File Exchange. 
%       [2] Şen, B., Peker, M., Çavuşoğlu, A., and Çelebi, F. V. (2014). A Comparative Study on Classification of Sleep Stage Based on EEG Signals Using Feature Selection and Classification Algorithms. J Med Syst 38, 18. doi: 10.1007/s10916-014-0018-0.
%       [3] Steven M. Pincus. Approximate Entropy as a Measure of System 
%           Complexity[J]. Proceedings of the National Academy of Sciences 
%           of the United States of America, 1991, 88(6):2297-2301. 
%           DOI: 10.1073/pnas.88.6.2297
%       [4] https://zhuanlan.zhihu.com/p/494761890  （code）
%       [5] https://zhuanlan.zhihu.com/p/574732876  （theory）

% 计算时间相对较长 Computationally time-consuming

function apen = feat_ApproximateEntropy(X, opts)

% 可选，通过参数tau对输入X降采样后再计算ApEn
% if nargin < 4, tau = 1; end
% if tau > 1, X = downsample(X, tau); end

X = X(:);  % 强制转化数据为列方向
if isfield(opts,'r')	% 指定
    r = opts.r; 
else	% 未指定 default
    r = 0.2;
end
r = r*std(X);
N = length(X);  % 信号长度，通常在100-5000范围内才能保证会有有效的统计特性和较小的误差
% dim:重构维数，序列长度为100-1000时，推荐dim=2，序列长度为1000-30000时，推荐dim=3
if (100<N) && (N<1000) 
    dim = 2;
else
    dim = 3;
end
phi = zeros(1,2);

for j = 1:2
    m = dim+j-1;        % dim通常设置为2，这样m=2和m=3; dim=3时，m=3和m=4
    C = zeros(1,N-m+1);   % 近似比例
    dataMat = zeros(m,N-m+1);   % 子序列集合，每个子序列的长度是m，共从data中提取N-m+1个子序列
    % setting up data matrix
    for i = 1:m
        dataMat(i,:) = X(i:N-m+i);
    end

    % counting similar patterns using distance calculation
    for i = 1:N-m+1
        tempMat = abs(dataMat - repmat(dataMat(:,i),1,N-m+1));
        boolMat = any( (tempMat > r),1);
        C(i) = sum(~boolMat)/(N-m+1);         % 近似比例
    end

    % summing over the counts
    phi(j) = sum(log(C))/(N-m+1);
end
apen = phi(1)-phi(2);
end