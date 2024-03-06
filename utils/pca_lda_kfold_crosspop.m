function [perf,Post_all] = pca_lda_kfold_crosspop(XX,YY,XX2,YY2,param)

ncd = length(unique(YY{1}));
Post_all = zeros(length(XX),ncd,ncd);
for t = 1 : length(XX)
    group=[YY{t}  ; YY2{t}];
    data_sub = [XX{t} ; XX2{t}] ;
    
  %  [coeff,score,~,~,explained,~] = pca(data_sub,'NumComponents',param.nComp);
  %  data_sub = score;

    %- TRAINING ON JUICE 1 and TESTING ON JUICE 2
    trainingData = data_sub(group<200,:);
    sampleData = data_sub(group>200,:);

    %- try classification.. failed when too many units have non-zeros FR
    class = [];
    try [class,err,posterior,logp,coeff] = classify(sampleData,trainingData,group(group<200), 'diaglinear');
    end

    if ~isempty(class)
        class = class+100;
        testi = group(group>200);
        cd=unique(group(group>200));
        for tr = 1 : length(class)
            Post_all(t,cd==class(tr),cd==testi(tr)) = Post_all(t,cd==class(tr),cd==testi(tr))+1;
        end
        perf(t,1) = mean(abs(class-group(group>200))==0);
    else
        perf(t,1) = NaN;
    end

    %- TRAINING ON JUICE 2 and TESTING ON JUICE 1
    trainingData = data_sub(group>200,:);
    sampleData = data_sub(group<200,:);

    %- try classification.. failed when too many units have non-zeros FR
    class = [];
    try [class,err,posterior,logp,coeff] = classify(sampleData,trainingData,group(group>200), 'diaglinear');
    end

    if ~isempty(class)
        class = class-100;
        testi = group(group<200);
        cd=unique(group(group<200));
        for tr = 1 : length(class)
            Post_all(t,cd==class(tr),cd==testi(tr)) = Post_all(t,cd==class(tr),cd==testi(tr))+1;
        end
        perf(t,2) = mean(abs(class-group(group<200))==0);
    else
        perf(t,2) = NaN;
    end


end

perf = nanmean(perf,2);



