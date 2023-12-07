%% 非线(神经动力学)性特征 | nonlinear features 
%% permutation Lempel-Ziv complexity, PLZC   （P219）
% 相比普通的LZC_mean、LZC_median、LZC_midpoint、LZC_kmeans，PLZC具有计算简单、鲁棒性/稳定性强、计算复杂度低等优点
% LZC的值越大，代表序列越复杂，越接近随机状态  
% X: single channel EEG signal (either a row vector or a column vector)
% opts.pattern：将信号表示为由不超过m!种模式组成的序列，如3,4

%% Reference
%       [1] Bai, Y., Li, X., Liang, Z. (2019). Nonlinear Neural Dynamics. In: Hu, L., Zhang, Z. (eds) EEG Signal Processing and Feature Extraction. Springer, Singapore. https://doi.org/10.1007/978-981-13-9113-2_11


function plzc = feat_PLempelZivComplexity(x, opts)
        
    if isfield(opts,'pattern')	% 指定
        m = opts.pattern; 
    else	% 未指定 default
        m = 4;
    end
    x = x(:)';  % 行向量
    datap=PX(x,m,1);    % 使用排序算法将脑电信号转换成有限数字序列，通过该步骤信号被表示为不超过m!种代表排序模式的符号
    
    m_factorial = factorial(m);     % m!
    c = 1;                                                                   
    S = datap(1);Q = [];
    for i=2:length(datap)
        Q = strcat(Q,datap(i));
        SQ = strcat(S,Q);
        SQv = SQ(1:length(SQ)-1);
        if isempty(findstr(SQv,Q))      % 如果Q不是SQv中的子串，说明Q是新出现的模式，c++    
            S = SQ;
            Q = [];
            c = c+1;    
        end
    end
    b = length(datap)*log(m_factorial)/log(length(datap));   
    plzc = c/b;
end

function OPi = PX(S,ord,t)

    ly = length(S);
    permlist = perms(1:ord);                    %列出1到ord的所有可能的排序
    c(1:length(permlist))=0;

    OPi=zeros(1,ly-t*(ord-1));                  

     for j=1:ly-t*(ord-1)
         [a,iv]=sort(S(j:t:j+t*(ord-1)));       %对S中的ord个数升序排列返回iv是a中各个数据排列的标号
         for jj=1:length(permlist)
             if (abs(permlist(jj,:)-iv))==0     %比较排序模式得到数据的排序模式
                 OPi(j)=jj;                     %OPi中是数据对应的排序模式
             end
         end
     end  
end