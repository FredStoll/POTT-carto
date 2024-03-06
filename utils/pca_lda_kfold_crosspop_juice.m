function [perf,Post_all] = pca_lda_kfold_crosspop_juice(XX,YY,XX2,YY2,param)

ncd = length(unique(YY{1}));
Post_all = zeros(length(XX),ncd,ncd);
for t = 1 : length(XX)
    group=[YY{t}  ; YY2{t}];
    data_sub = [XX{t} ; XX2{t}] ;
    
%   [coeff,score,~,~,explained,~] = pca(data_sub,'NumComponents',param.nComp);
%   data_sub = score;

   XX{1} = data_sub(ismember(group,unique(YY{1})),:);
   XX2{1} = data_sub(ismember(group,unique(YY2{1})),:);

    %- train on 1st test on second
    %- try classification.. failed when too many units have non-zeros FR
    class = [];
    try [class,err,posterior,logp,coeff] = classify(XX2{1},XX{1},YY{1}, 'diaglinear');
    end

    if ~isempty(class)
        class = floor(class/100);
        testi = floor(YY2{1}/100);
        cd=unique(group(group>200));
        for tr = 1 : length(class)
            Post_all(t,cd==class(tr),cd==testi(tr)) = Post_all(t,cd==class(tr),cd==testi(tr))+1;
        end
        perf(t,1) = mean(abs(class-testi)==0);
    else
        perf(t,1) = NaN;
    end

    %- opposite classifier
    %- try classification.. failed when too many units have non-zeros FR
    class = [];
    try [class,err,posterior,logp,coeff] = classify(XX{1},XX2{1},YY2{1}, 'diaglinear');
    end

    if ~isempty(class)
        class = floor(class/100);
        testi = floor(YY{1}/100);
        cd=unique(group(group>200));
        for tr = 1 : length(class)
            Post_all(t,cd==class(tr),cd==testi(tr)) = Post_all(t,cd==class(tr),cd==testi(tr))+1;
        end
        perf(t,2) = mean(abs(class-testi)==0);
    else
        perf(t,2) = NaN;
    end





end


%perf = nanmean(perf,2);



