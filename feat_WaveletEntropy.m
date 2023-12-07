%% 非线性（动力学）特征：基于时频的熵 | nonlinear features: time-frequency-based entropy
%% 小波熵 | Wavelet Entropy (WE) (P226)
% 1.计算每一段信号尺度的小波能量
% 2.将各小波能量除以总能量，以获得各尺度j的相对小波能量
% 3.WE即不同尺度间分布的熵
% 当脑电信号表现出更强的规律性时，即神经元活动出现同步的现象，小波熵值将降低
% WE的值取决于小波基函数、分解层数n和数据长度N
% 其中，小波基函数在小波分析中最为重要，由于缺乏固定的准则，实际应用中通常很难直接选择合适的小波基函数，通常都是根据实验结果来确定该函数
% X: single channel EEG signal (either a row vector or a column vector)
% opts.fs: sampling frequency 信号采样率
% opts.wavefunc: wavelet basis function, default = 'db4'


%% Reference
%       [1] Bai, Y., Li, X., Liang, Z. (2019). Nonlinear Neural Dynamics. In: Hu, L., Zhang, Z. (eds) EEG Signal Processing and Feature Extraction. Springer, Singapore. https://doi.org/10.1007/978-981-13-9113-2_11


function wentropy = feat_WaveletEntropy(X, opts)
    
    fs = opts.fs;
    if isfield(opts,'wavefunc')	% 指定
        wavefunc = opts.wavefunc; 
    else	% 未指定 default
        wavefunc = 'db4';
    end
    
    Nlayer=round(log2(fs))-3;   % 分解层数，fs=200时Nlayer=5
    wentropy=0;
    E=waveletdecom_cwq(X,Nlayer,wavefunc);     % 这里默认选择了db4小波基函数
    P=E/sum(E);
    P = P(find(P~=0));
    for j=1:size(P,2)
        wentropy=wentropy-P(1,j).*log(P(1,j));    %小波熵Swt=-sum(Pj*logPj)
    end

end

function [ E ] = waveletdecom_cwq( x,n,wpname )
[C,L]=wavedec(x,n,wpname);        %对数据进行小波包分解
for k=1:n
      %wpcoef(wpt1,[n,i-1])是求第n层第i个节点的系数
%       disp('每个节点的能量E(i)');
      SRC(k,:)=wrcoef('a',C,L,'db4',k);%尺度
      SRD(k,:)=wrcoef('d',C,L,'db4',k);%细节系数

  E(1,k)=norm(SRD(k,:))*norm(SRD(k,:));
      %求第i个节点的范数平方，其实也就是平方和
end
E(1,n+1)=norm(SRC(n,:))*norm(SRC(n,:));
%  disp('小波包分解总能量E_total');
E_total=sum(E);  %求总能量
y=E_total;
% disp('以下是每个节点的概率P');
% for i=1:n+1
%    p(i)= E(1,i)/E_total    %求每个节点的概率
% end

end
