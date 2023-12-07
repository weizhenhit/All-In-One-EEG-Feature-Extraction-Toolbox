%% 包络变异系数（函数输入通常是imfx而非X） | Variation Coefficient of the Envelope (CVE)(based on imfs)
% 基于解析表示的瞬时振幅，反映包络离散程度的绝对值
% 1. 调用hilbert()函数对信号进行Hilbert变换得到信号的解析表示
% 2. 基于信号的解析表示，计算瞬时振幅（即信号包络）和瞬时相位
% 3. 基于信号包络计算包络变异系数CVE = std(envelope)/mean(envelope)
% imfx: IMF component of single channel EEG signal (either a row vector or
% a column vector)  当然直接基于单通道EEG信号计算其谱统计特征也可以
% opts.fs: sampling frequency 信号采样率
% opts.window_size: 滑动窗口大小，default=1s

%% Reference
%       [1] Díaz, J., Bassi, A., Coolen, A., Vivaldi, E. A., and Letelier, J.-C. (2018). Envelope analysis links oscillatory and arrhythmic EEG activities to two types of neuronal synchronization. NeuroImage 172, 575–585. doi: 10.1016/j.neuroimage.2018.01.063.）


function CVE = feat_CoefficientVariationEnvelope(imfx, opts)
    fs = opts.fs;
    if isfield(opts,'window_size')	% 指定
        window_size = opts.window_size; 
    else	% 未指定 default
        window_size = 1;
    end
    
    % X = hilbert(Xr) computes the so-called discrete-time analytic signal
    % X = Xr + i*Xi such that Xi is the Hilbert transform of real vector Xr.
    imfx_ana = hilbert(imfx);   % 对信号进行Hilbert变换得到信号的解析表示,length=N
    envelope = abs(imfx_ana);   % 信号的瞬时振幅instantaneous amplitude，即信号包络，length=N
    % 将包络信号分割成窗口
    % window_size = 1;    % 窗口长度：1s
    num_windows = floor(length(envelope) / (fs * window_size));
    envelope_windows = reshape(envelope(1:num_windows*fs*window_size), fs*window_size, num_windows);
    % 计算每个窗口的变异系数 CEV=std(envelope)/mean(envelope)
    cv = std(envelope_windows, 0, 1) ./ mean(envelope_windows, 1);      
    % 取所有窗口的变异系数的平均值作为包络变异系数
    CVE = mean(cv);
end