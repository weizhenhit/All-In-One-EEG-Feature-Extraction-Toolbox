% 根据白化法或目标函数法构造出的空间滤波器矩阵CSPMatrix，提取CSP特征

function features = func_extractCSPFeatures_seg(EEG, CSPMatrix, nbFilterPairs,seg_len)
% Input
%       EEG: 训练集或测试集
%           EEG.x: (times, channels, trials)
%           EEG.y: trials*1
%       CSPMetrix: 学习到的CSP空间滤波器矩阵, cell, {1,nsegments},每个都是channels*channels
%       nbFilterPairs: m
%       seg_len: 每个segment的长度，points
% Return
%       features: CSP特征矩阵,(trials, 2m, nsegments)

nbTrials = size(EEG.x,3);  
nsegments = length(CSPMatrix);
features = zeros(nbTrials, 2*nbFilterPairs, nsegments);    

for t=1:nbTrials   
    trial = EEG.x(:,:,t);   % (times*channels)
    for seg=1:nsegments
        % 每个segment的CSPMatrix只能作用于每个segment上，而不能作用于整个trial
        M = CSPMatrix{1, seg};  % channels*channels
        Filter = M([1:nbFilterPairs (end-nbFilterPairs+1):end],:);  % 2m * channels
        % 1:nbFilterPairs是前m=nbFilterPairs行，(end-nbFilterPairs+1):end是后m=nbFilterPairs行
        
        segment = trial((seg-1)*seg_len+1:seg*seg_len, :);      % (times*channels)
        projectedSeg = Filter * segment';    % Filter:(2m*channels), EEG segment:(times*channels)，projectedTrial:(2m*times)
        % 投影信号的对数方差作为生成特征
        variances = var(projectedSeg,0,2);    % length = 2m

        for f=1:length(variances)
            % features(t,f) = variances(f)/sum(variances); 
            % features(t,f) = log(variances(f)/sum(variances)); 
            % 下面一种的acc更高一些，所以使用下面的实现方式
            features(t,f,seg) = log(variances(f));      % (nbTrials, 2*nbFilterPairs, nsegments)
        end
    end
end
