% 根据白化法或目标函数法构造出的空间滤波器矩阵CSPMatrix，提取CSP特征

function features = func_extractCSPFeatures(EEG, CSPMatrix, nbFilterPairs)
% Input
%       EEG: 训练集或测试集
%           EEG.x: (times, channels, trials)
%           EEG.y: trials*1
%       CSPMetrix: 学习到的CSP空间滤波器矩阵, 每个都是channels*channels
%       nbFilterPairs: m
% Return
%       features: CSP特征矩阵,(trials, 2m)

nbTrials = size(EEG.x,3);  
features = zeros(nbTrials, 2*nbFilterPairs);    
Filter = CSPMatrix([1:nbFilterPairs (end-nbFilterPairs+1):end],:);  
% 1:nbFilterPairs是前m=nbFilterPairs行，(end-nbFilterPairs+1):end是后m=nbFilterPairs行，所以得到的Filter是2m * channels的矩阵

for t=1:nbTrials    
    projectedTrial = Filter * EEG.x(:,:,t)';    % Filter是2m*channels, EEG trial矩阵X是times*channels，得到的projectedTrial是2m*times
    % 投影信号的对数方差作为生成特征
    variances = var(projectedTrial,0,2);    % length = 2m

    for f=1:length(variances)
        % features(t,f) = variances(f)/sum(variances); 
        % features(t,f) = log(variances(f)/sum(variances)); 
        % 下面一种的acc更高一些，所以使用下面的实现方式
        features(t,f) = log(variances(f));      % features: CSP特征矩阵,(trials, 2m)
    end
%     features(t,end) = EEGSignals.y(t);                            %不要拼接标签   最后一列为类标签
end
