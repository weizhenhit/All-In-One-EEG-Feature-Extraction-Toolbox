%% 时频特征 | time-frequency domain
%% 经验模态分解 | empirical mode decomposition (EMD)
% 使用 EMD 将复杂的信号分解并简化为执行希尔伯特谱分析所需的有限个固有模态函数/内涵模态分量（intrinsic mode
% functions，IMFs）和一个残余分量residual，每一个IMF代表了原信号不同频率段的振荡变化，反映信号的局部特征，而最后的残余分量则反映信号中的缓慢变化量。
% IMFs函数的特性：
% ①在整个数据段内，极值点的个数和过零点的个数必须相等或相差最多不能超过一个
%       即图线要反复跨越x轴，而不能在某次穿过零点后出现多个极点
% ②在任意时刻，由局部极大值点形成的上包络线和由局部极小值点形成的下包络线的平均值为零
%       即上、下包络线相对于时间轴局部对称。
% EMD works like an adaptive high pass filter. 
% It shifts out the fastest changing component first and as the level of IMF increases, 
% the oscillation of IMF becomes smoother. Each component is band-limited, 
% which can reflect the characteristic of instantaneous frequency.
% X: single channel EEG signal (either a row vector or a column vector)
% opts.fs: sampling frequency 信号采样率
% opts.imf_level: only use the first 'imf_level' imfs for later feature extraction

%% Reference
%       [1] Mert, A., and Akan, A. (2018). Emotion recognition from EEG signals by using multivariate empirical mode decomposition. Pattern Anal Appl 21, 81–89. doi: 10.1007/s10044-016-0567-6.
%       [2] Riaz, F., Hassan, A., Rehman, S., Niazi, I. K., and Dremstrup, K. (2016). EMD-Based Temporal and Spectral Features for the Classification of EEG Signals Using Supervised Learning. IEEE Transactions on Neural Systems and Rehabilitation Engineering 24, 28–35. doi: 10.1109/TNSRE.2015.2441835.
%       [3] https://zhuanlan.zhihu.com/p/40005057
%       [4] https://blog.csdn.net/qq_40061206/article/details/120664537

%% emd()函数示例
% load('sinusoidalSignalExampleData.mat','X','fs')  %载入数据
% t = (0:length(X)-1)/fs;
% plot(t,X)                       %绘制原始信号图:由不同振幅和频率的正弦信号叠加得到
% xlabel('Time(s)') 
% emd(X,'Interpolation','pchip')  %emd分解
% % 从title中可以看到，一种有9个IMF分量，而图中只显示了其中的IMF1~IMF3，如果要显示其他分量，在图片的空白处点击右键
% [imf, residual, info] = emd(x);
% % 函数返回值主要包括imf\residual\info三个，imf即各模态分量值；residual维残差值；info中包括了该次分解的一些重要信息，如imf数量、各分量的过零点数、各分量的极值数等

function imf_features = feat_EmpiricalModeDecomposition(X, opts)
    fs = opts.fs;
    if isfield(opts,'imf_level')	% 指定
        imf_level = opts.imf_level;     % Maximum number of IMFs extracted, one of the decomposition stop criteria
    else	% 未指定 default
        imf_level = 3;      % 默认最大提取前3个IMFs
    end
        
    %% EMD
    [imf, residual, info] = emd(X, 'Interpolation','pchip','MaxNumIMF',imf_level ,'Display',1);     % 对于非平稳信号使用pchip插值
    % imf: [N, imf_level]
    
    %% IMFs分量的选择
    % 可以只选择2阶IMF（即IMF2）进行下一步的特征提取
    % 一方面，2阶IMF频率基本集中在0-60Hz，是正常脑电能量集中的频带；
    % 另一方面，通过1阶IMF初步滤掉了高频噪声，2阶IMF的信号质量有所改善
    % 但是考虑到后续的特征提取和分类步骤，特征数开始尽可能多，之后再通过特征选择（例如选择通道数）减少特征维度
    % 默认imf_level = 3，表示只提取前三个IMFs
    
    %% 特征提取
    % 通常会选取若干阶IMFs，从这些IMFs中提取这些信号的时域和频域特征，作为脑电信号的时频域特征
    % 1. 时域：直接基于N阶IMFs时域信号，计算其各种时域特征
            % 统计特征（平均值、方差、偏度、峰度、一阶差分平均值）、
            % Hjorth特征（Hjorth Activity、Hjorth Mobility、Hjorth Complexity）
            % 非线性特征（分形维数、mean energy）
            % 共10 * imf_level * channels (太多了，imf_level得取小一些N=3，特征/通道选择+降维)
    % 2. 频域：首先对N阶IMFs时域信号进行功率谱密度估计得pxx，再提取基于pxx的频谱特征
            % 频谱质心、频谱方差、频谱偏度、频谱峰度
            % 共4 * imf_level * channels（特征/通道选择+降维）
    % 3. 解析表示：对N阶IMFs时域信号进行Hilbert变换得到IMFs的解析表示，提取其解析表示的统计特征等
            % 包络变异系数（算频域特征，基于解析表示的瞬时振幅，反映包络离散程度的绝对值：Díaz, J., Bassi, A., Coolen, A., Vivaldi, E. A., and Letelier, J.-C. (2018). Envelope analysis links oscillatory and arrhythmic EEG activities to two types of neuronal synchronization. NeuroImage 172, 575–585. doi: 10.1016/j.neuroimage.2018.01.063.）
            % 共1 * imf_level * channels（特征/通道选择+降维）
            
    % todo:IMFAR_features = feat_IMFsAnalyticRepresentation(X, opts)
    % IMFs的解析表示：对N阶IMFs时域信号进行Hilbert变换得到IMFs的解析表示，提取其解析表示的统计特征等
    
    imf_features = [];
    for imf_num = 1:imf_level
        imfx = imf(:, imf_num);
        %% 1. 时域特征
        fmean = feat_Mean(imfx);
        fvar = feat_Variance(imfx);
        fske = feat_Skewness(imfx);
        fkur = feat_Kurtosis(imfx);
        f1d = feat_FirstDifference(imfx);
        fhjac = feat_HjorthActivity(imfx);
        fhjmo = feat_HjorthMobility(imfx);
        fhjco = feat_HjorthComplexity(imfx);
        ffrdi = feat_FractalDimension(imfx);
        fmeen = feat_MeanEnergy(imfx);
        %% 2. 频域特征
        opts_.fs = fs;
        [fcentroid, fvariance, fskewness, fkurtosis] = feat_SpectrumStatistical(imfx, opts_);   
        
        %% 3. 解析表示
        % X = hilbert(Xr) computes the so-called discrete-time analytic signal
        % X = Xr + i*Xi such that Xi is the Hilbert transform of real vector Xr.
        % imfx_ana = hilbert(imfx);   % 对信号进行Hilbert变换得到信号的解析表示,length=N
        % envelope = abs(imfx_ana);   % 信号的瞬时振幅instantaneous amplitude，即信号包络，length=N
        % 包络变异系数CVE
        opts_.window_size = 1;  % 滑动窗口长度=1s
        fcve = feat_CoefficientVariationEnvelope(imfx, opts_);
        
        feature_imfx = [fmean fvar fske fkur f1d fhjac fhjmo fhjco ffrdi fmeen fcentroid fvariance fskewness fkurtosis fcve];
        imf_features = [imf_features feature_imfx];
    end
end