%% 非线性（动力学）特征：基于时域的熵 | nonlinear features: time-domain-based entropy 
%% 排列熵 | Permutation Entropy 
% 表征时间序列数据的复杂性（complexity）和规律性（regularity），通过计算序列中不同排列模式的概率分布来衡量序列的不确定性，值越大，时间序列越复杂、越不规律
% X: single channel EEG signal (either a row vector or a column vector)
% opts.tau: 延迟时间 delay time
% See also:
%       feat_ShannonEntropy, feat_SampleEntropy, feat_ApproximateEntropy, feat_FuzzyEntropy

%% Reference
%       [1] Bandt C，Pompe B. Permutation entropy:a natural complexity measure for time series[J]. Physical Review Letters,2002,88(17):174102.
%       [2] https://blog.csdn.net/weixin_45317919/article/details/109254213
%       [3] https://blog.csdn.net/weixin_45317919/article/details/124516421

function peEn = feat_PermutationEntropy(X,opts)

if isfield(opts,'tau')	% 指定
    t = opts.tau; 
else	% 未指定 default
    t = 1;
end

X = X(:)';      % 强制转化数据为行向量
N = length(X);
% dim:重构维数，序列长度为100-1000时，推荐dim=2，序列长度为1000-30000时，推荐dim=3
if (100<N) && (N<1000) 
    m = 2;
else
    m = 3;
end
permlist = perms(1:m);
[h,~]=size(permlist);
c(1:length(permlist))=0;

 for j=1:N-t*(m-1)
     [~,iv]=sort(X(j:t:j+t*(m-1)));
     for jj=1:h
         if (abs(permlist(jj,:)-iv))==0
             c(jj) = c(jj) + 1 ;
         end
     end
 end
hist = c;
c=c(c~=0);
p = c/sum(c);
peEn = -sum(p .* log(p));
% 归一化
peEn=peEn/log(factorial(m));
end
