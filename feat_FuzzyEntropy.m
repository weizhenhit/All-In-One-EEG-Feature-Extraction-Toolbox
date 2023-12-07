%% 非线性（动力学）特征：基于时域的熵 | nonlinear features: time-domain-based entropy 
%% 模糊熵 | Fuzzy entropy （based on Sample Entropy）
% 表征时间序列数据的复杂性（complexity）和规律性（regularity），值越大，时间序列越复杂、越不规律
% X: single channel EEG signal (either a row vector or a column vector)
% opts.r：相似容限阈值，通常r在0.1std(X)和0.25std(X)之间时，En有比较合理的统计特性，默认0.2std(X)
% opts.n:决定了相似容限边界的梯度，n越大则梯度越大，n在计算向量间相似性中起着权重的作用，为了捕获尽量多的细节信息，提出者建议计算时取较小的整数，如2,3
% See also:
%       feat_ShannonEntropy, feat_SampleEntropy, feat_ApproximateEntropy, feat_PermutationEntropy

%% Reference
%       [1] Chen Weiting,Wang Zhizhong,XieHongbo,et al. Characterization of surfaceEMG signal based on fuzzy entropy. IEEE Transactions on Neural Systems and Rehabilitation Engineering. 2007,15(2):266-272.
%       [2] https://zhuanlan.zhihu.com/p/519633544  （code）
%       [3] https://zhuanlan.zhihu.com/p/574793525  （theory）
%       [4] https://blog.csdn.net/weixin_45317919/article/details/110008733


function FuzEn = feat_FuzzyEntropy(X, opts)
    if isfield(opts,'r')	% 指定
        r = opts.r; 
    else	% 未指定 default
        r = 0.2;
    end
    if isfield(opts,'n')	% 指定
        n = opts.n; 
    else	% 未指定 default
        n = 2;
    end
    
    X = X(:)';      % 强制转化数据为行向量
    N = length(X);
    % dim:重构维数，序列长度为100-1000时，推荐dim=2，序列长度为1000-30000时，推荐dim=3
    if (100<N) && (N<1000) 
        dim = 2;
    else
        dim = 3;
    end
    r = r * std(X);
    phi = zeros(1, 2);

    for m = dim:dim+1
       count = zeros(N-m+1, 1);
       dataMat = zeros(N-m+1, m);
       
       % 设置数据矩阵，构造成dim维的矢量
       for i = 1:N-m+1
          dataMat(i, :) = X(1, i:i+m-1)-mean(X(1, i:i+m-1));   % (1)序列定义 
       end
       % 利用距离计算相似模式数
       for j = 1:N-m+1
          % 计算切比雪夫距离，不包括自匹配情况
          tempmat = repmat(dataMat(j,:), N-m+1, 1);
          dist = max(abs(dataMat - tempmat), [], 2);    % (2)遍历所有X(i)和X(j)的组合，并求出端点坐标差的绝对值的最大值
          D = exp(-(dist.^n)/r);        % (3)模糊隶属度
          count(j) = (sum(D)-1)/(N-m-1);      % （4）分子中-1是为了减掉i=j的情况
       end
       phi(m-dim+1) = sum(count)/(N-m);   % (4)
    end
    FuzEn = log(phi(1)/phi(2));     % (5)
end