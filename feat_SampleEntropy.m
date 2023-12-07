%% 非线性（动力学）特征：基于时域的熵 | nonlinear features: time-domain-based entropy 
%% 样本熵 | Sample Entropy （similar to Approximate Entropy）
% 表征时间序列数据的复杂性（complexity）和规律性（regularity），值越大，时间序列越复杂、越不规律
% X: single channel EEG signal (either a row vector or a column vector)
% opts.r：相似容限阈值，通常r在0.1std(X)和0.25std(X)之间时，En有比较合理的统计特性，默认0.2std(X)
% See also:
%       feat_ShannonEntropy, feat_ApproximateEntropy, feat_PermutationEntropy, feat_FuzzyEntropy

%% Reference
%       [1] Kijoon Lee (2023). Sample Entropy (https://www.mathworks.com/matlabcentral/fileexchange/35784-sample-entropy), MATLAB Central File Exchange.
%       [2] Richman JS, Moorman JR. 
%           Physiological time-series analysis using approximate entropy and 
%           sample entropy[J]. 
%           Am J Physiol Heart Circ Physiol, 2000, 278(6):H2039-49. 
%           DOI: 10.1152/ajpheart.2000.278.6.H2039.
%       [3] https://zhuanlan.zhihu.com/p/519042619  （code）
%       [4] https://zhuanlan.zhihu.com/p/574764496  （theory）

% 计算时间相对较长 Computationally time-consuming

function saen = feat_SampleEntropy(X, opts)

%   SampleEntropy is conceptually similar to approximate entropy (ApEn), but has
%   following differences:
%       1) SampEn does not count self-matching. The possible trouble of
%       having log(0) is avoided by taking logarithm at the latest step.
%       2) SampEn does not depend on the datasize as much as ApEn does. The
%       comparison is shown in the graph that is uploaded.

% if nargin < 4, tau = 1; end
% if tau > 1, X = downsample(X, tau); end

X = X(:);   % 强制转化数据为列方向
if isfield(opts,'r')	% 指定
    r = opts.r; 
else	% 未指定 default
    r = 0.2;
end
r = r*std(X);
N = length(X);
% dim:重构维数，序列长度为100-1000时，推荐dim=2，序列长度为1000-30000时，推荐dim=3
if (100<N) && (N<1000) 
    dim = 2;
else
    dim = 3;
end
phi = zeros(1,2);

for m = dim:dim+1
    Bi = zeros(1,N-m+1);
    dataMat = zeros(m,N-m+1);

    % setting up data matrix
    for i = 1:m
        dataMat(i,:) = X(i:N-m+i);
    end

    % counting similar patterns using distance calculation
    for j = 1:N-m+1
        % calculate Chebyshev distance, excluding self-matching case
        dist = max(abs(dataMat - repmat(dataMat(:,j),1,N-m+1)));
        % calculate Heaviside function of the distance
        % User can change it to any other function
        % for modified sample entropy (mSampEn) calculation
        D = (dist <= r);
        % excluding self-matching case
        Bi(j) = (sum(D)-1)/(N-m);
    end

    % summing over the counts
    phi(m-dim+1) = sum(Bi)/(N-m+1);

end

saen = -log(phi(2)/phi(1));

end