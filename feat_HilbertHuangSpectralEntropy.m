%% 非线性（动力学）特征：基于时频的熵 | nonlinear features: time-frequency-based entropy
%% 希尔伯特-黄谱熵 | Hilbert-Huang Spectral Entropy (HHSE) (P227)
% 1. 对于非平稳时域信号x(t)，经EMD将其分解为一系列IMFs（1,...,N）
% 2. 对IMFs分量进行Hilbert变换，得到IMF的解析表示
% 3. 基于IMFs的解析表示，计算归一化后的Hilbert-Huang spectral
% 4. 将香农熵算法应用到Hilbert-Huang spectral中，即得到Hilbert-Huang spectral entropy
% 注意：输入是单通道EEG时域信号，最终输出的HHSE特征只有一个数值，而非每个IMF信号对应一个HHSE特征
% X: single channel EEG signal (either a row vector or a column vector)
% opts.fs: sampling frequency 信号采样率
% opts.f_bin: frequend range [f_low f_high], e.g.[1 45]

% called functions：emd2，extr，io，hhspectrum，toimage，instfreq

%% Reference
%       [1] Bai, Y., Li, X., Liang, Z. (2019). Nonlinear Neural Dynamics. In: Hu, L., Zhang, Z. (eds) EEG Signal Processing and Feature Extraction. Springer, Singapore. https://doi.org/10.1007/978-981-13-9113-2_11


function [En, P, HHTP] = feat_HilbertHuangSpectralEntropy(X, opts)

% return:  En - entropy
%          P -  marginal spectrum 边际谱：指在某一特定范围内的光谱，通常位于光谱的边缘部分，可能包含一些不太明显或不太重要的信息
% input:   X - time series
%          f_bin - frequend band [f_low, f_high]
%          fs - sampling frequency
%          df - frequency resolution 默认=fs/1000,即nfft=1000

f_bin = opts.f_bin;
fs = opts.fs;
df = fs/1000;

X = reshape(X, 1, length(X));
X = [X X X];    % 重复了3次  似乎与im = im(:,N/3:2*N/3-1);有关

N=length(X);    % disp(['N=', num2str(N)]);

%% EMD
[imf,ort,nbits] = func_emd2(X);

k = size(imf);      % disp('size of imf='); disp(k);
rx=zeros(1,N);

%  imf = imf(1:end-1,:);

%% Hilbert spectrum
[A,f,tt] = func_hhspectrum(imf); % disp(size(A));  % （10，3454）  3454=1152*3-2

[im,t]=func_toimage(A,f,tt,fs/2/df);     % transforms a spectrum made of 1D functions (e.g., output of "spectreh") in an 2D image
% show Hilbert spectrum
% disp_hhs(im,t);

im = im(:,N/3:2*N/3-1);     % 恢复到长度=N

HHTP = flipud(im);
% marginal spectrum
P=sum(flipud(im),2)/N;      % Flip翻转 array in up/down direction (upside down).

Pf = P(fix(f_bin(1)/df)+1:fix(f_bin(2)/df));    % fix(X) rounds the elements of X to the nearest integer towards zero. 如fix(-1.9)=-1, fix(1.9)=1
En = en_log(Pf);

function en = en_log(P)

P=P/sum(P);
P = P(find(P~=0));
en = 0;
for i = 1:length(P)
    en = en - P(i)*log(P(i));
end

