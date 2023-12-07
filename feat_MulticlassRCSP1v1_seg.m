%% 共空间模式 | common spatial patterns (CSP)  不可在main中调用，仅可单独使用
%% one vs one Multiclass or two-classes  目标函数法
% how：N个类别中任意两个做CSP，共得到CN2=N*(N-1)/2个投影矩阵，对每个投影矩阵取最佳的2m个方向（前后各m个），特征矩阵为2m*N*(N-1)/2=m*N*(N-1)维
% 优点：每次计算CSP投影矩阵时，2个类别都是单纯的类别，更容易找到两类差别最大的投影方向
% 缺点：N类需要计算CN2个投影矩阵，N很大时，计算量会大大增加
% 注意：在计算CSP投影矩阵时，并非使用整个trial，而是将trial划分出若干segments，
%   对每个segments提取CSP特征，最终只选择分类效果最好的（或者特定任务下某个特别选择的）segment或者所有segments的CSP特征作为该trial的CSP特征
%   此时，每个trial的特征矩阵表示为(m*N*(N-1), segments)
%   以(Lee et al., 2023)为例，每个trial(2s, 2500Hz, 5000points)分成了16个非重叠的segments(0.125s,312.5points)，每个segments都提取CSP特征，使用1vR策略，m=4
%   最终得到的特征矩阵是104 features * 16 segments，其中104 features = 13 classes * 8 csp
% Inputs
%       EEG_train: 训练集各类别的所有trials：
%           EEG_train.x: EEG，times*channels*trials
%           EEG_train.y: 训练集所有样本的类别标签，trials*1
%       EEG_test: 测试集各类别的所有trials：
%           EEG_test.x: EEG，times*channels*trials
%           EEG_test.y: 测试集所有样本的类别标签，trials*1
%       opts.fs: sampling frequency
%       opts.m: number of CSP filters, dimensions of CSP features are 2m
%       opts.seg_len: length of each segments，sec
% Return
%       CSPMatrix: cell，学习到的CSP滤波器矩阵，共{nbClasses*(nbClasses-1)/2, nsegments}个，每个都是(channels * channels)
%       feature_train: 
%           feature_train.data: 训练集对应的CSP特征矩阵，(trials, m*nbClasses*(nbClasses-1), nsegments)
%           feature_train.y: 即EEG_train.y
%       feature_test: 
%           feature_test.data: 测试集对应的CSP特征矩阵，(trials, m*nbClasses*(nbClasses-1), nsegments)
%           feature_test.y: 即EEG_test.y
% called function：func_extractCSPFeatures_seg
% See also
%       feat_MulticlassCSP1v1_seg, feat_MulticlassCSP1vR_seg,
%       feat_MulticlassRCSP1vR_seg

%% Reference
%       [1] https://blog.csdn.net/qq_40166660/article/details/115218031
%       [2] Naeem, M., Brunner, C., Leeb, R., Graimann, B., and Pfurtscheller, G. (2006). Seperability of four-class motor imagery data using independent components analysis. J. Neural Eng. 3, 208. doi: 10.1088/1741-2560/3/3/003.
%       [3] 刘广权,黄淦,朱向阳.共空域模式方法在多类别分类中的应用[J].中国生物医学工程学报,2009,28(06):935-938.
%       [4] https://blog.csdn.net/MissXy_/article/details/81264953
%       [5] https://blog.csdn.net/oh__NO/article/details/84310982


% 基于最优化问题构造空间滤波器矩阵 RSCP(目标函数法)
function [feature_train, feature_test, CSPMatrix] = feat_MulticlassRCSP1v1_seg(EEG_train, EEG_test, opts)
    fs = opts.fs;
    m = opts.m;
    if isfield(opts,'seg_len')	% 指定
        seg_len = opts.seg_len;     % 每个segment的长度，s
    else	% 未指定 default
        seg_len = 1;
    end
    seg_len = floor(seg_len * fs);  % 每个trial的长度，points
    
    ntimes = size(EEG_train.x, 1);  % 每个trial的times数
    nbChannels = size(EEG_train.x,2);        
    nbTrials = size(EEG_train.x,3);          
    classLabels = unique(EEG_train.y);       % 类标签，如[1,2,3,4]
    nbClasses = length(classLabels);          % 类别数
    disp([num2str(nbClasses) ' classes csp!']);
    nsegments = floor(ntimes/seg_len);      % segments数
    
    %计算每个trial的每个segment标准化后的协方差矩阵
    trialCov = zeros(nbChannels,nbChannels,nbTrials,nsegments);
    for t = 1:nbTrials
        E = EEG_train.x(:,:,t)';    % x: times*channels*trials  E: channels*times
        for seg = 1:nsegments
            E_seg = E(:, (seg-1)*seg_len+1:seg*seg_len);   
            EE = E_seg * E_seg';    % channels*channels
            trialCov(:,:,t, seg) = EE ./ trace(EE);
        end
    end
    clear E; clear E_seg; clear EE;
    
    % 每次选取N个类别中的任意2类
    nbMetrix = nbClasses * (nbClasses-1) / 2;   % 共CN2种情况
    CSPMatrix = cell(nbMetrix, nsegments);
    feature_train_data = zeros(nbTrials, 2*m*nbMetrix, nsegments);  % 每个trial的特征矩阵为（2*m*nbMetrix, nsegments）即每个样本可看做2维而非1维
    feature_test_data = zeros(size(EEG_test.x,3), 2*m*nbMetrix, nsegments);     % 每个trial的label
    feature_train.y = EEG_train.y;      % label
    feature_test.y = EEG_test.y;
    
    count = 0;
    for c1=1:nbClasses-1
        for c2 = c1+1:nbClasses
            count = count + 1;
            for seg = 1:nsegments

                covMatrices1 = mean(trialCov(:,:,EEG_train.y == classLabels(c1),seg),3);    % c1类
                covMatrices1 = covMatrices1 ./ sum(EEG_train.y == classLabels(c1));     % 存疑，是否应该再除一遍
                covMatrices2 = mean(trialCov(:,:,EEG_train.y == classLabels(c2),seg),3);    % c2类
                covMatrices2 = covMatrices2 ./ sum(EEG_train.y == classLabels(c2));

                %% 以下开始使用RCSP，基于Lagrangian乘子法，将CSPMatrix的构造问题转换为最优化问题
                % 计算需要分解的矩阵M  a\b 等价与 inv(a) * b
                M = covMatrices2 \ covMatrices1;

                %矩阵M的特征值分解
                [U D] = eig(M);
                eigenvalues = diag(D);
                [eigenvalues egIndex] = sort(eigenvalues, 'descend');
                U = U(:,egIndex);
                % 得到空间滤波器矩阵
                CSPMatrix{count, seg} = U';      % nbChannels * nbChannels
                % disp(size(CSPMatrix{count, seg}));
            end
            %% 提取CSP特征  训练集和测试集都使用训练集中得到的CSP空间滤波器矩阵
            % feature_train_2m和feature_test_2m的shape都是(trials,2m,nsegments)
            feature_train_2m = func_extractCSPFeatures_seg(EEG_train,CSPMatrix(count, :),m,seg_len);     %提取训练集特征       
            feature_test_2m = func_extractCSPFeatures_seg(EEG_test,CSPMatrix(count, :),m,seg_len);       %提取测试集特征
            % disp(size(feature_train_2m));
            % disp(size(feature_test_2m));
            
            % feature_train_data和feature_test_data的shape都是(trials, 2m*CN2, nsegments)
            feature_train_data(:, (count-1)*2*m+1:count*2*m, :) = feature_train_2m;     
            feature_test_data(:, (count-1)*2*m+1:count*2*m, :) = feature_test_2m;
        end
    end
    feature_train.data = feature_train_data;    % 每个trial的特征矩阵为（2*m*nbMetrix, nsegments）即每个样本可看做2维而非1维
    feature_test.data = feature_test_data;
end
