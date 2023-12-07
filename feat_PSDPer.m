%% 功率谱密度 | power spectral density
%% 非参数方法-周期图法 | Periodogram
% 直接调用bandpower(X, fs, [f_low, f_high])完成功率谱估计和相对average band power特征提取
% 功率谱使用默认的modified periodogram函数估计
% X: single channel EEG signal (either a row vector or a column vector)
% opts.fs: sampling frequency 信号采样率
% opts.f_low: freqrange — Frequency range for band power computation
% opts.f_high: freqrange — Frequency range for band power computation
% return: 5频带和2Hz均分频带划分得到的（5+bandwidth/2）维特征
% See also:
%       feat_PSDWelch, feat_PSDArburg

%% Reference
%       [1] Zhang, Z. (2019). Spectral and Time-Frequency Analysis. In: Hu, L., Zhang, Z. (eds) EEG Signal Processing and Feature Extraction. Springer, Singapore. https://doi.org/10.1007/978-981-13-9113-2_6
%       [2] https://github.com/zhangzg78/eegbook
%       [3] https://blog.csdn.net/weixin_44425788/article/details/127036076
%       [4] https://blog.csdn.net/frostime/article/details/106967703


function PSD = feat_PSDPer(X, opts)
    X = detrend(X);     % EEG频谱分析之前通常需要detrend去趋势
    fs = opts.fs;
    f_low = opts.f_low;
    f_high = opts.f_high;
    
    N = length(X);   % signal的长度
    nfft = 2^nextpow2(N); % the number of FFT points   FFT点数一般取大于信号点数N的最小的二次幂
    
    band_power_all = bandpower(X, fs, [f_low f_high]);
    
    % 传统5频带
    f_low_fre = [f_low, 4, 8, 12, 30];
    f_high_fre = [4, 8, 12, 30, f_high];
    average_band_power_5 = zeros(1, length(f_low_fre));
    for num = 1:length(f_low_fre)
        average_band_power_5(1, num) = bandpower(X, fs, [f_low_fre(num) f_high_fre(num)]);  % (1, bands=5) 
        average_band_power_5(1, num) = average_band_power_5(1, num)/band_power_all;
    end
    % disp(average_band_power_5);
    
    % 2Hz均分频带划分
    bands = (f_high-f_low)/2;
    average_band_power_2 = zeros(1, bands);
    for num = 1:bands
        average_band_power_2(1, num) = bandpower(X, fs, [(num-1)*2+1 num*2+1]);  % (1, bands=22) 
        average_band_power_2(1, num) = average_band_power_2(1, num)/band_power_all; 
    end
    % disp(average_band_power_2);
    
    % 其他特征：将pxx看做一个随机过程，从中提取其统计值（如均值、方差、熵等）作为频谱特征
    PSD = [average_band_power_5 average_band_power_2];
end