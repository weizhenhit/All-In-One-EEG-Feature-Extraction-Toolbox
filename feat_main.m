%% 脑电特征提取工具箱 | All-In-One EEG Feature Extraction Toolbox 
% An all-in-one EEG feature extraction toobox, including statistical features, 
%   Hjorth parameters, entropy, nonlinear features, power spectral density (PSD), 
%   differential entropy (DE), empirical mode decomposition (EMD), common spatial patterns (CSP), 
%   microstate analysis and so on. (The list of features will continue to update...)
%% 1.时域特征 | time domain （channels*1）
%   1.1.统计特征 | statistical features
%       (1) Mean    'mean'    @feat_Mean
%       (2) Standard Deviation      'sd'    @feat_StandardDeviation
%       (3) Variance   'var'   @feat_Variance
%       (4) Median     'md'    @feat_Median
%       (5) Maximum    'max'   @feat_Maximum
%       (6) Minimum    'min'   @feat_Minimum
%       (7) First Difference/Mean Curve Length    '1d'    @feat_FirstDifference
%       (8) Normalized First Difference     'n1d'   @feat_NormalizedFirstDifference
%       (9) Second Difference   '2d'    @feat_SecondDifference
%       (10) Normalized Second Difference    'n2d'   @feat_NormalizedSecondDifference
%       (11) Kurtosis     'kurt'  @feat_Kurtosis
%       (12) Skewness     'skew'  @feat_Skewness
%       (13) Zero-crossing Rate     'c0r'   @feat_ZeroCrossingRate
%       (14) Log Root Sum Of Sequential Variation     'lrssv'    @feat_LRSSV
%   1.2.Hjorth parameters
%       (15) Hjorth Activity     'ha'    @feat_HjorthActivity
%       (16) Hjorth Mobility     'hm'    @feat_HjorthMobility
%       (17) Hjorth Complexity   'hc'    @feat_HjorthComplexity
%   1.3.非线性特征：基于时域的熵 | nonlinear features: time-domain-based entropy
%       (18) Shannon Entropy    'sh'    @feat_ShannonEntropy
%       (19) Sample Entropy         'se'    @feat_SampleEntropy
%       (20) Approximate Entropy    'ae'    @feat_ApproximateEntropy
%       (21) Permutation Entropy    'pe'    @feat_PermutationEntropy
%       (22) Fuzzy Entropy      'fe'    @feat_FuzzyEntropy
%       (23) Renyi Entropy      're'    @feat_RenyiEntropy
%       (24) Tsallis Entropy    'te'    @feat_TsallisEntropy
%       (25) Log Energy Entropy     'le'    @feat_LogEnergyEntropy
%   1.4.非线性特征 | nonlinear features
%       (26) Mean Energy        'me'        @feat_MeanEnergy
%       (27) Mean Teager Energy     'mte'   @feat_MeanTeagerEnergy
%       (28) Fractal Dimension      'fd'    @feat_FractalDimension
%       (29) Permutation Lempel-Ziv Complexity      'plzc'       @feat_PLempelZivComplexity 
%       (30) Recurrence Quantification Analysis    'rqa'    @feat_RecurrenceQA
%       (31) Hurst Exponent        'he'    @feat_HurstExponent
%% 2.频域特征 | frequency domain
%   2.1.功率谱估计+power spectral density特征提取（特定频带内的相对average band power：2Hz均分频带法、5频带法，同时返回两种方式计算得到的PSD特征）  channels*bands
%       (32) Periodogram      'psd_per'   @feat_PSDPer  （直接调用bandpower(X, fs, [f_low, f_high])提取特征,默认使用modified periodogram方法估计功率谱密度pxx）
%       (33) Welch      'psd_welch'         @feat_PSDWelch
%       (34) Auto Regressive Model (burg)   'psd_arburg'    @feat_PSDArburg  （分别调用pwelch和pburg估计功率谱密度pxx，然后调用bandpower(pxx, f, [f_low, f_high], 'psd')提取特征）
%       非参数方法(32 33)：信号包含大量噪声、信号长度足够长(可以理解为划分的窗口比较多,平均之后更加平滑)时的结果更准确
%       参数方法(34)：若模型准确，即使对于较短的信号也可以得到可靠的频谱估计，对噪声和模型参数敏感
%   2.2.微分熵 | differential entropy (DE)
%       (35) Differential Entropy   'de'    @feat_DifferentialEntropy
%   2.3.基于IMFs的谱统计特征（函数输入通常是imfx而非X） | spectral statistical features (based on imfs)
%       (36) Spectral centroid      'specs'     @feat_SpectrumStatistical
%       (37) Spectral variance      'specs'     @feat_SpectrumStatistical
%       (38) Spectral skewness      'specs'     @feat_SpectrumStatistical
%       (39) Spectral kurtosis      'specs'     @feat_SpectrumStatistical
%   2.4.包络变异系数 | Variation Coefficient of the Envelope (CVE)
%       (40) Variation Coefficient of the Envelope   'cve'    @feat_CoefficientVariationEnvelope
%% 3.时频特征 | time-frequency domain
%   3.1.经验模态分解 | empirical mode decomposition (EMD)
%       (41) Empirical Mode Decomposition   'emd'   @feat_EmpiricalModeDecomposition
%   3.2.非线性（动力学）特征：基于时频的熵 | nonlinear features: time-frequency-based entropy
%       (42) Hilbert-Huang Spectral Entropy     'hhse'   @feat_HilbertHuangSpectralEntropy
%       (43) Wavelet Entropy   'we'     @feat_WaveletEntropy
%% 4.空域特征 | spatial domain
%   4.1.共空间模式 | common spatial patterns (CSP)  不可在main中调用，仅可单独使用
%       (44) Common Spatial Patterns (白化法 one vs. one)   'csp1v1'   @feat_MulticlassCSP1v1_seg   default
%       (45) Common Spatial Patterns (白化法 one vs. rest)   'csp1vR'   @feat_MulticlassCSP1vR_seg
%       (46) Common Spatial Patterns (目标函数法 one vs. one)   'rcsp1v1'  @feat_MulticlassRCSP1v1_seg
%       (47) Common Spatial Patterns (目标函数法 one vs. rest)   'rcsp1vR'  @feat_MulticlassRCSP1vR_seg
%   4.2.微状态分析 | microstate analysis   不可在main中调用，仅可单独使用
%       (48)Microstate Analysis      

%% Inputs
% X: single channel EEG signal (either a row vector or a column vector)
% opts: necessary or optional parameters for each feature extraction function (refer to the 'Parameter’column of [List of available features])

%% Outputs
% feature: extracted feature vector or feature matrix (the dimension of `feature` is listed in‘Dimensions' column of [List of available features])

%% Examples: demo.m
% You can call and test each feature extraction function separately.
% You can also use the 'feat_main' function to extract multiple features at once. The idea of feat_main.m comes from https://github.com/JingweiToo/EEG-Feature-Extraction-Toolbox
% The usage of multiclass Common Spatial Patterns function 'feat_MulticlassCSP1v1_seg' is different from others. Details can be found in feat_MulticlassCSP1v1_seg.m

%% Reference
%       [1] Zhang, Z. (2019). Spectral and Time-Frequency Analysis. In: Hu, L., Zhang, Z. (eds) EEG Signal Processing and Feature Extraction. Springer, Singapore. https://doi.org/10.1007/978-981-13-9113-2_6
%       [2] Jia, H. (2019). Microstate Analysis. In: Hu, L., Zhang, Z. (eds) EEG Signal Processing and Feature Extraction. Springer, Singapore. https://doi.org/10.1007/978-981-13-9113-2_8
%       [3] Bai, Y., Li, X., Liang, Z. (2019). Nonlinear Neural Dynamics. In: Hu, L., Zhang, Z. (eds) EEG Signal Processing and Feature Extraction. Springer, Singapore. https://doi.org/10.1007/978-981-13-9113-2_11
%       [4] Jenke, R., Peer, A., and Buss, M. (2014). Feature Extraction and Selection for Emotion Recognition from EEG. IEEE Transactions on Affective Computing 5, 327–339. doi: 10.1109/TAFFC.2014.2339834.
%       [5] https://github.com/JingweiToo/EEG-Feature-Extraction-Toolbox
%       [6] https://github.com/zhangzg78/eegbook
%       [7] https://blog.csdn.net/weixin_44425788/category_12026253.html


function feature = feat_main(abbr, X, opts)
switch abbr
    case 'mean'   ; fun = @feat_Mean;
    case 'sd'     ; fun = @feat_StandardDeviation;
    case 'var'    ; fun = @feat_Variance;
    case 'md'     ; fun = @feat_Median;
    case 'max'    ; fun = @feat_Maximum;
    case 'min'    ; fun = @feat_Minimum;
    case '1d'     ; fun = @feat_FirstDifference;
    case 'n1d'    ; fun = @feat_NormalizedFirstDifference;
    case '2d'     ; fun = @feat_SecondDifference;
    case 'n2d'    ; fun = @feat_NormalizedSecondDifference;
    case 'kurt'   ; fun = @feat_Kurtosis;
    case 'skew'   ; fun = @feat_Skewness;
    case 'c0r'    ; fun = @feat_ZeroCrossingRate;
    case 'lrssv'  ; fun = @feat_LRSSV;   
    case 'ha'     ; fun = @feat_HjorthActivity;
    case 'hm'     ; fun = @feat_HjorthMobility;
    case 'hc'     ; fun = @feat_HjorthComplexity;
    case 'sh'     ; fun = @feat_ShannonEntropy;
    case 'se'     ; fun = @feat_SampleEntropy;               opts.r = 0.2; 
    case 'ae'     ; fun = @feat_ApproximateEntropy;          opts.r = 0.2;
    case 'pe'     ; fun = @feat_PermutationEntropy;          opts.tau = 1;
    case 'fe'     ; fun = @feat_FuzzyEntropy;                opts.r = 0.2;   opts.n = 2;     
    case 're'     ; fun = @feat_RenyiEntropy;                opts.alpha = 2;
    case 'te'     ; fun = @feat_TsallisEntropy;              opts.alpha = 2;
    case 'le'     ; fun = @feat_LogEnergyEntropy;
    case 'me'     ; fun = @feat_MeanEnergy;
    case 'mte'    ; fun = @feat_MeanTeagerEnergy;    
    case 'fd'     ; fun = @feat_FractalDimension;
    case 'plzc'   ; fun = @feat_PLempelZivComplexity;       opts.pattern = 4;     % 将信号表示为由不超过m!种模式组成的序列，如3,4
    case 'rqa'    ; fun = @feat_RecurrenceQA;                opts.dimension = 3;     opts.delay = 1;     opts.thr = 0.5;     opts.index = 0;
    case 'he'     ; fun = @feat_HurstExponent;               % opts.fs = ;    opts.start_frequency = ;   opts.step_frequency = ;     opts.stepnumber = ;
    case 'psd_per'; fun = @feat_PSDPer;                      % opts.fs = ;    opts.f_low = ;    opts.f_high = ;        
    case 'psd_welch'; fun = @feat_PSDWelch;                  % opts.fs = ;    opts.f_low = ;    opts.f_high = ;        
    case 'psd_arburg'; fun = @feat_PSDArburg;                opts.order = 20;    % opts.fs = ;     opts.f_low = ;    opts.f_high = ;    
    case 'de'     ; fun = @feat_DifferentialEntropy;        
    case 'specs'  ; fun = @feat_SpectrumStatistical;         % opts.fs = ;      % 包含频谱质心、频谱方差、频谱偏度、频谱峰值、频道总能量5种频谱统计特征
    case 'cve'    ; fun = @feat_CoefficientVariationEnvelope;   opts.window_size = 1;   % opts.fs = ;         % 滑动窗口大小=1s
    case 'emd'    ; fun = @feat_EmpiricalModeDecomposition;  % opts.fs = ;      opts.imf_level = 3;   % 为避免特征过多，这里仅使用前3个IMFs
    case 'hhse'   ; fun = @feat_HilbertHuangSpectralEntropy; % opts.fs = ;      opts.f_bin = [1 45];    
    case 'we'     ; fun = @feat_WaveletEntropy;              opts.wavefunc = 'db4';      % opts.fs = ;
    case 'csp1v1' ; fun = @feat_MulticlassCSP1v1_seg;       % opts.m = ;   opts.fs = ;     opts.seg_len = 1;   % seg_len单位为s，默认segments长度为1s
    case 'csp1vR' ; fun = @feat_MulticlassCSP1vR_seg;       % opts.m = ;   opts.fs = ;     opts.seg_len = 1;
    case 'rcsp1v1'; fun = @feat_MulticlassRCSP1v1_seg;      % opts.m = ;   opts.fs = ;     opts.seg_len = 1;
    case 'rcsp1vR'; fun = @feat_MulticlassRCSP1vR_seg;      % opts.m = ;   opts.fs = ;     opts.seg_len = 1;
end
if nargin < 3
    opts = [];
end
feature = fun(X, opts);
end

