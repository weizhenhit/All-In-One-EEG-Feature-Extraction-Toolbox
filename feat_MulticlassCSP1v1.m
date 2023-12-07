%% 共空间模式 | common spatial patterns (CSP)  不可在main中调用，仅可单独使用
%% one vs one Multiclass or two-classes 白化法
% how：N个类别中任意两个做CSP，共得到CN2=N*(N-1)/2个投影矩阵，对每个投影矩阵取最佳的2m个方向（前后各m个），特征矩阵为2m*N*(N-1)/2=m*N*(N-1)维
% 优点：每次计算CSP投影矩阵时，2个类别都是单纯的类别，更容易找到两类差别最大的投影方向
% 缺点：N类需要计算CN2个投影矩阵，N很大时，计算量会大大增加
% Inputs
%       EEG_train: 训练集各类别的所有trials：
%           EEG_train.x: EEG，times*channels*trials
%           EEG_train.y: 训练集所有样本的类别标签，trials*1
%       EEG_test: 测试集各类别的所有trials：
%           EEG_test.x: EEG，times*channels*trials
%           EEG_test.y: 测试集所有样本的类别标签，trials*1
%       m: number of CSP filters, dimensions of CSP features are 2m
% Return
%       CSPMatrix: cell，学习到的CSP滤波器矩阵，共nbClasses * (nbClasses-1) / 2个，每个都是(channels * channels)
%       feature_train: 训练集对应的CSP特征矩阵，(trials, 1+m*nbClasses*(nbClasses-1))
%       feature_test:  测试集对应的CSP特征矩阵，(trials, 1+m*nbClasses*(nbClasses-1))
% called function：func_extractCSPFeatures
% See also
%       feat_MulticlassCSP1vR, feat_MulticlassRCSP1v1, feat_MulticlassRCSP1vR

%% Reference
%       [1] https://blog.csdn.net/qq_40166660/article/details/115218031
%       [2] Naeem, M., Brunner, C., Leeb, R., Graimann, B., and Pfurtscheller, G. (2006). Seperability of four-class motor imagery data using independent components analysis. J. Neural Eng. 3, 208. doi: 10.1088/1741-2560/3/3/003.
%       [3] 刘广权,黄淦,朱向阳.共空域模式方法在多类别分类中的应用[J].中国生物医学工程学报,2009,28(06):935-938.
%       [4] https://blog.csdn.net/MissXy_/article/details/81264953
%       [5] https://blog.csdn.net/oh__NO/article/details/84310982

% 基于传统特征向量和白化矩阵构造空间滤波器矩阵（白化法）
function [feature_train, feature_test, CSPMatrix] = feat_MulticlassCSP1v1(EEG_train, EEG_test, m)
    nbChannels = size(EEG_train.x,2);        
    nbTrials = size(EEG_train.x,3);          
    classLabels = unique(EEG_train.y);       % 类标签，如[1,2,3,4]
    nbClasses = length(classLabels);          % 类别数
    disp([num2str(nbClasses) ' classes csp!']);
    
%     covMatrices = cell(nbClasses,1);  %每一类的协方差矩阵
    %计算每个trial标准化后的协方差矩阵
    trialCov = zeros(nbChannels,nbChannels,nbTrials);
    for t=1:nbTrials
        E = EEG_train.x(:,:,t)';    % times*channels*trials
        EE = E * E';
        trialCov(:,:,t) = EE ./ trace(EE);
    end
    clear E;
    clear EE;
    
    % 每次选取N个类别中的任意2类
    nbMetrix = nbClasses * (nbClasses-1) / 2;   % 共CN2种情况
    CSPMatrix = cell(nbMetrix, 1);
    feature_train = zeros(nbTrials, 1+2*m*nbMetrix);
    feature_test = zeros(size(EEG_test.x,3), 1+2*m*nbMetrix);
    feature_train(:, 1) = EEG_train.y;      % 返回的特征矩阵的第一列是label列
    feature_test(:, 1) = EEG_test.y;
    
    count = 0;
    for c1=1:nbClasses-1
        for c2 = c1+1:nbClasses
            count = count + 1;
            % 计算2类的平均协方差矩阵
    %         for c=1:nbClasses
    %             covMatrices{c} = mean(trialCov(:,:,EEGSignals.y == classLabels(c)),3);
    %         end
            covMatrices1 = mean(trialCov(:,:,EEG_train.y == classLabels(c1)),3);    % c1类
            covMatrices2 = mean(trialCov(:,:,EEG_train.y == classLabels(c2)),3);    % c2类

            % 总的协方差矩阵 （原理中的混合空间协方差矩阵R）
    %         covTotal = covMatrices{1} + covMatrices{2};
            covTotal = covMatrices1 + covMatrices2;

            % 总协方差矩阵的白化转换(eig求特征值和特征向量)
            [Ut Dt] = eig(covTotal);     % 注意：特征值初始是升序排列的
            eigenvalues = diag(Dt);
            [eigenvalues egIndex] = sort(eigenvalues, 'descend');
            Ut = Ut(:,egIndex);
            P = diag(sqrt(1./eigenvalues)) * Ut';     % 白化特征矩阵

            % 通过P进行第一类协方差矩阵的转换   
            % 可证明transformedCov1与transformedCov2的特征向量矩阵是相同的，且可证明transformedCov1的最大特征值所对应的特征向量使transformedCov2有最小的特征值，反之亦然
            % 所以后续只使用transformedCov1即可，二者是一样的
            transformedCov1 =  P * covMatrices1 * P';

            % 协方差矩阵的EVD
            % 对transformedCov1分解得特征向量和对应特征值
            [U1 D1] = eig(transformedCov1);
            eigenvalues = diag(D1);
            [eigenvalues egIndex] = sort(eigenvalues, 'descend');       % 把transformedCov1的特征值按降序排序（对应transformedCov1的特征值按升序排序），在最大的和最小的特征值中各选择m个，这2m个对应的特征向量构成了U1
            U1 = U1(:, egIndex);
            % 得到空间滤波器矩阵
            CSPMatrix{count} = U1' * P;    
            % disp(size(CSPMatrix{count}));     % nbChannels * nbChannels

            %% 提取CSP特征  训练集和测试集都使用训练集中得到的CSP空间滤波器矩阵
            % feature_train和feature_test的shape都是(trials,2m)，2m代表每个样本的CSP特征维数是2m=4，与通道数无关
            feature_train_2m = func_extractCSPFeatures(EEG_train,CSPMatrix{count},m);     %提取训练集特征       
            feature_test_2m = func_extractCSPFeatures(EEG_test,CSPMatrix{count},m);       %提取测试集特征
            % disp(size(feature_train_2m));
            % disp(size(feature_test_2m));

            feature_train(:, 1+(count-1)*2*m+1:1+count*2*m) = feature_train_2m;
            feature_test(:, 1+(count-1)*2*m+1:1+count*2*m) = feature_test_2m;
        end
    end
end

