%% 共空间模式 | common spatial patterns (CSP)  不可在main中调用，仅可单独使用
%% one vs rest Multiclass or two-classes  目标函数法
% how：N个类别得到N个投影矩阵，对每个投影矩阵取最佳的2m个方向（前后各m个），特征矩阵为2m*N维
% 优点：N类需要计算N个投影矩阵，即使N很大，计算量也不会很大
% 缺点：由于rest类别是混合类别，内部特征并不统一，因此很难找出最优的投影方向，即结果可能并不是最优的
% Inputs
%       EEG_train: 训练集各类别的所有trials：
%           EEG_train.x: EEG，times*channels*trials
%           EEG_train.y: 训练集所有样本的类别标签，trials*1
%       EEG_test: 测试集各类别的所有trials：
%           EEG_test.x: EEG，times*channels*trials
%           EEG_test.y: 测试集所有样本的类别标签，trials*1
%       m: number of CSP filters, dimensions of CSP features are 2m
% Return
%       CSPMatrix: cell，学习到的CSP滤波器矩阵，共nbClasses个，每个都是(channels * channels)
%       feature_train: 训练集对应的CSP特征矩阵，(trials, 2m*nbClasses)
%       feature_test:  测试集对应的CSP特征矩阵，(trials, 2m*nbClasses)
% called function：func_extractCSPFeatures
% See also
%       feat_MulticlassCSP1v1, feat_MulticlassCSP1vR, feat_MulticlassRCSP1v1

%% Reference
%       [1] https://blog.csdn.net/qq_40166660/article/details/115218031
%       [2] Naeem, M., Brunner, C., Leeb, R., Graimann, B., and Pfurtscheller, G. (2006). Seperability of four-class motor imagery data using independent components analysis. J. Neural Eng. 3, 208. doi: 10.1088/1741-2560/3/3/003.
%       [3] 刘广权,黄淦,朱向阳.共空域模式方法在多类别分类中的应用[J].中国生物医学工程学报,2009,28(06):935-938.
%       [4] https://blog.csdn.net/MissXy_/article/details/81264953
%       [5] https://blog.csdn.net/oh__NO/article/details/84310982


% 基于最优化问题构造空间滤波器矩阵 RSCP(目标函数法)
function [feature_train, feature_test, CSPMatrix] = feat_MulticlassRCSP1vR(EEG_train, EEG_test, m)
    nbChannels = size(EEG_train.x,2);        
    nbTrials = size(EEG_train.x,3);          
    classLabels = unique(EEG_train.y);       % 类标签，如[1,2,3,4]
    nbClasses = length(classLabels);          % 类别数
    disp([num2str(nbClasses) ' classes csp!']);
    
%     covMatrices = cell(nbClasses,1);  %每一类的协方差矩阵
    %计算每个trial标准化后的协方差矩阵
    trialCov = zeros(nbChannels,nbChannels,nbTrials);
    for t=1:nbTrials
        E = EEG_train.x(:,:,t)';
        EE = E * E';
        trialCov(:,:,t) = EE ./ trace(EE);
    end
    clear E;
    clear EE;
    
    % 只有两类：当前类与剩余类
    if nbClasses == 2
        CSPMatrix = cell(1, 1);
        feature_train = zeros(nbTrials, 1+2*m);
        feature_test = zeros(size(EEG_test.x,3), 1+2*m); 
    else
        CSPMatrix = cell(nbClasses, 1);
        feature_train = zeros(nbTrials, 1+2*m*nbClasses);
        feature_test = zeros(size(EEG_test.x,3), 1+2*m*nbClasses);
    end
    feature_train(:, 1) = EEG_train.y;      % 返回的特征矩阵的第一列是label列
    feature_test(:, 1) = EEG_test.y;
    
    for c=1:nbClasses
        
        % 计算2类的平均协方差矩阵
%         for c=1:nbClasses
%             covMatrices{c} = mean(trialCov(:,:,EEGSignals.y == classLabels(c)),3);
%         end
        covMatrices1 = mean(trialCov(:,:,EEG_train.y == classLabels(c)),3);    % 当前类
        covMatrices1 = covMatrices1 ./ sum(EEG_train.y == classLabels(c));
        covMatrices2 = mean(trialCov(:,:,EEG_train.y ~= classLabels(c)),3);    % 剩余类
        covMatrices2 = covMatrices2 ./ sum(EEG_train.y ~= classLabels(c));
        
        %% 以下开始使用RCSP，基于Lagrangian乘子法，将CSPMatrix的构造问题转换为最优化问题
        % 计算需要分解的矩阵M  a\b 等价与 inv(a) * b
        M = covMatrices2 \ covMatrices1;

        %矩阵M的特征值分解
        [U D] = eig(M);
        eigenvalues = diag(D);
        [eigenvalues egIndex] = sort(eigenvalues, 'descend');
        U = U(:,egIndex);
        % 得到空间滤波器矩阵
        CSPMatrix{c} = U';      % nbChannels * nbChannels
        disp(size(CSPMatrix{c}));
        
        %% 提取CSP特征  训练集和测试集都使用训练集中得到的CSP空间滤波器矩阵
        % feature_train和feature_test的shape都是(trials,2m)，2m代表每个样本的CSP特征维数是2m=4，与通道数无关
        feature_train_2m = func_extractCSPFeatures(EEG_train,CSPMatrix{c},m);     %提取训练集特征       
        feature_test_2m = func_extractCSPFeatures(EEG_test,CSPMatrix{c},m);       %提取测试集特征
        disp(size(feature_train_2m));
        disp(size(feature_test_2m));
        
        feature_train(:, 1+(c-1)*2*m+1:1+c*2*m) = feature_train_2m;
        feature_test(:, 1+(c-1)*2*m+1:1+c*2*m) = feature_test_2m;
        
        if nbClasses == 2
            % 如果只有2个类别，没必要算两遍
            return
        end
    end
end
