%% 基于IMFs的谱统计特征（函数输入通常是imfx而非X） | spectral statistical features (based on imfs)
% Spectral centroid 频谱质心: 信号频率分布和能量分布的重要信息，是频率成分的重心，是在一定频率范围内通过能量加权平均的频率
% Spectral variance 频谱方差: 反映频域能量变化的波动程度
% Spectral skewness 频谱偏度: 表征频谱在质心附近的不对称性
% Spectral kurtosis 频谱峰度: 表征频谱在质心附近的平坦度
% sum band power 频带能量
% EMD分解得到的IMFs通常具有窄带特点，每个IMFs代表了原信号不同频率段的振荡变化
% 默认使用welch方法[pxx,f] = pwelch(imfx,fs,fs/2,nfft,fs)进行功率谱估计
% imfx: IMF component of single channel EEG signal (either a row vector or
% a column vector)  当然直接基于单通道EEG信号计算其谱统计特征也可以
% opts.fs: sampling frequency 信号采样率
% return: 5种基于IMFx的谱统计特征

%% Reference
%       [1] Riaz, F., Hassan, A., Rehman, S., Niazi, I. K., and Dremstrup, K. (2016). EMD-Based Temporal and Spectral Features for the Classification of EEG Signals Using Supervised Learning. IEEE Transactions on Neural Systems and Rehabilitation Engineering 24, 28–35. doi: 10.1109/TNSRE.2015.2441835.

function [centroid, variance, skewness, kurtosis, band_power] = feat_SpectrumStatistical(imfx, opts)
    % centroid: 频谱质心
    % variance: 频谱方差
    % skewness: 频谱偏度
    % kurtosis: 频谱峰度
    % band_power: sum band power
    
    fs = opts.fs;
    N = length(imfx);   % signal的长度
    nfft = 2^nextpow2(N); % the number of FFT points   FFT点数一般取大于信号点数N的最小的二次幂
    imfx = detrend(imfx);     % EEG频谱分析之前通常需要detrend去趋势
    [pxx, f] = pwelch(imfx,fs,fs/2,nfft,fs);       % [pxx,f] = pwelch(x,window,noverlap,nfft,fs)
    
    % 计算频谱质心
    % The spectral centroid is the first order moment (一阶矩) of the frequency based on the energy distribution, 
    % which can reflect the characteristic of the fundamental harmonic (基谐波) of EEG. 
    % In addition, the smaller value of spectral centroid indicates that more energy is concentrated on the low frequency range
    centroid = sum(f .* pxx) / sum(pxx);

    % 计算频谱方差
    variance = sum((f - centroid).^2 .* pxx) / sum(pxx);

    % 计算频谱偏度
    skewness = sum(((f - centroid) / sqrt(variance)).^3 .* pxx) / sum(pxx);
    
    % 计算频谱峰度
    kurtosis = sum(((f - centroid) / sqrt(variance)).^4 .* pxx) / sum(pxx);
    
    % 计算band power  未指定f_low和f_high表示计算所有f频带上的band power
    % 因为考虑到每个IMFs似乎对应一段特定的频率成分，所以用5频带法似乎不太合适，直接用全频带比较合适
    band_power = bandpower(pxx, f, 'psd');              % p = bandpower(pxx,f,'psd')
end
