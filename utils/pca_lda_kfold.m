function [perf,Post_all] = pca_lda_kfold(XX,YY,param)

ncd = length(unique(YY{1}));
Post_all = zeros(length(XX),ncd,ncd);
for t = 1 : length(XX)
    group=YY{t};
    data_sub = XX{t};
    
    [coeff,score,~,~,explained,~] = pca(data_sub,'NumComponents',param.nComp);
    data_sub = score;

    cv = cvpartition(group,'KFold',10); % 0.20'
    numFolds = cv.NumTestSets;
    for i = 1:numFolds
        sampleInds = cv.test(i);
        trainInds = cv.training(i);

        trainingData = data_sub(trainInds,:);
        sampleData = data_sub(sampleInds,:);
        
        %- try classification.. failed when too many units have non-zeros FR
        class = [];
        try [class,err,posterior,logp,coeff] = classify(sampleData,trainingData,group(trainInds), 'diaglinear');
        end

        if ~isempty(class)
            testi = group(sampleInds);
            cd=unique(group);
            for tr = 1 : length(class)
                Post_all(t,cd==class(tr),cd==testi(tr)) = Post_all(t,cd==class(tr),cd==testi(tr))+1;
            end
            mean(abs(class-group(sampleInds)));
            perf(t,i) = mean(abs(class-group(sampleInds))==0);
        else
            perf(t,i) = NaN;
        end

        
    end
end
perf = nanmean(perf,2);
