
function [ALL,param] = behav_model_proba(filename,perm,showfig)

%- load behav variable
load(filename)

%- some utils
take = @(arg) strcmp(TrialType_header,arg);
norm = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;
when = @(arg,cd) TrialType(strcmp(TrialType_header,arg),:)==cd ;

%- permutations
if perm
    
    pav = when('task',1); %- leave the PAV trial alone
    ntr4perm = sum(~pav);
    tr_order = randperm(ntr4perm); %- perm the others (ins + dyn)
    
    param2perm = {'TrialType' 't_FPon' 't_FPfix_first' 't_respfix' 't_resp'};
    for pp = 1 : length(param2perm)
        eval(['dumm = ' param2perm{pp} '(:,~pav);']);
        eval([param2perm{pp} '=[' param2perm{pp} '(:,pav),dumm(:,tr_order)];']);
    end
    
    %- reboot the function to send the new TrialType matrix!
    when = @(arg,cd) TrialType(strcmp(TrialType_header,arg),:)==cd ;
end


%- trials to consider
diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
restr = when('task',2) & when('brks',0) & diff_juice; %- take only INS trials (task == 2), completed (brks == 0) and diff juice

%- variable to predict: ChooseJ1 - 1 for J1 / 0 for J2
choice = abs(TrialType(take('I_chosenjuice'),restr)-2);

if length(choice) > 100  %- trash sessions with less than 100 trials diff juice
    
    %% Step 1 - GLM on Choices in different juice trials
    %- extract Proba of both options
    probaJ1( when('I_juiceL',1) & restr ) = TrialType(take('I_probaL'),when('I_juiceL',1) & restr);
    probaJ1( when('I_juiceR',1) & restr ) = TrialType(take('I_probaR'),when('I_juiceR',1) & restr);
    probaJ1 = probaJ1(restr);
    probaJ2( when('I_juiceL',2) & restr ) = TrialType(take('I_probaL'),when('I_juiceL',2) & restr);
    probaJ2( when('I_juiceR',2) & restr ) = TrialType(take('I_probaR'),when('I_juiceR',2) & restr);
    probaJ2 = probaJ2(restr);
    
    pbJ1 = norm(probaJ1/100); %- normalize proba
    pbJ2 = norm(probaJ2/100);
    
    %- extract Reaction time
    rt = t_respfix-t_resp(1,:);
    rt = norm(log(rt(restr))); %- log then min-max normalization (-1 to 1)
    rt(~isfinite(rt)) = NaN;
    
    %- extract Initiation time
    ft = t_FPfix_first-t_FPon(1,:);
    ft = norm(log(ft(restr))); %- log then min-max normalization (-1 to 1)
    ft(~isfinite(ft)) = NaN;
    
    %- extract prevJuice (0 if chose J1 before / 1 if chose J2 before) - this is not side related
    %- prevJuice == which juice he 'could have get' on the last completed trial (not necessarily
    %- rewarded).
    %- WARNING!!! FOR NOW, INCLUDE Dynamic task trials -- If trial before is Dynamic task, which juice he chose there (Morbier) or was offered (Mimic)
    
    currjuice = NaN(1,length(TrialType));
    currjuice(when('task',1)) = TrialType(take('P_juice'),when('task',1))-1;
    currjuice(when('task',2)) = TrialType(take('I_chosenjuice'),when('task',2))-1;
    currjuice(when('task',5)) = TrialType(take('D_chosenjuice'),when('task',5))-1;
    currjuice(when('task',3)) = TrialType(take('D_chosenjuice'),when('task',3))-1;
    
    %- put the last trials' juice when brks at current (or multiple breaks in a row, take the last one)
    noJuice = find(isnan(currjuice));
    for tr = 1 : length(noJuice)
        if noJuice(tr)>1
            currjuice(noJuice(tr)) = currjuice(noJuice(tr)-1);
        end
    end
    
    prevJuice = [NaN currjuice(1:end-1)];
    prevJuice = prevJuice(restr);
    prevJuice(prevJuice==0)=-1; %- a way to standardize it (-1 vs 1)
    
    %- extract the last received juice
    received = when('rew',1);
    currjuice(~received)=NaN; %- remove all the juice info from unrewarded trials
    
    %- put the last trials' RECEIVED juice when no rew or brks (or multiple breaks in a row, take the last one)
    %- no information on how far was this last reward
    noRwd = find(isnan(currjuice));
    for tr = 1 : length(noRwd)
        if noRwd(tr)>1
            currjuice(noRwd(tr)) = currjuice(noRwd(tr)-1);
        end
    end
    
    prevJuiceReceived = [NaN currjuice(1:end-1)];
    prevJuiceReceived = prevJuiceReceived(restr);
    prevJuiceReceived(prevJuiceReceived==0)=-1; %- a way to standardize it (-1 vs 1)
    
    %- define model
    % modelspec = 'choice ~ 1 + probJ1*probJ2 + ft + prevJuice' ;
    modelspec = 'choice ~ 1 + prob + prevJuice' ;
    % modelspec = 'choice ~ 1 + probJ1*probJ2 + prevJuiceReceived' ;
    % modelspec = 'choice ~ 1 + probJ1*probJ2 + rt' ;
    
    %- create Table and set categorical variables
    T = table(choice',log(probaJ1./probaJ2)',ft',rt',prevJuice',prevJuiceReceived','VariableNames',{'choice','prob','ft','rt','prevJuice','prevJuiceReceived'});
    T.choice=categorical(T.choice);
    T.prevJuice=categorical(T.prevJuice);
    T.prevJuiceReceived=categorical(T.prevJuiceReceived);
    
    lastwarn('', '');
    mdl = fitglm(T,modelspec,'Distribution','binomial','Link','logit')
    

  %- check if glm converged
    [~, warnId] = lastwarn();
    if strcmp(warnId,'stats:glmfit:IterationLimit')
        converge=0;
    else
        converge=1;
    end

    %- predict choices (for that, prevJuice is set as 0 (J1) and RT is median(rt)
    clear newf newc
    probas = [10:1:90];
    allpb = log(repmat(probas,size(probas')) ./ repmat(probas,size(probas'))');

        tab = table(allpb(:),nanmedian(ft)*ones(length(probas)*length(probas),1),nanmedian(rt)*ones(length(probas)*length(probas),1),-ones(length(probas)*length(probas),1),-ones(length(probas)*length(probas),1),'VariableNames',{'prob','ft','rt','prevJuice','prevJuiceReceived'});
        tab.prevJuice = categorical(tab.prevJuice);
        tab.prevJuiceReceived = categorical(tab.prevJuiceReceived);
    
    [newf , newc] = predict(mdl,tab);

    %% Extract param of interest for each session!
    ALL.name = filename;
    ALL.mdl = mdl;
    ALL.converge = converge;
    ALL.newf = newf;
    ALL.newc = newc;
    ALL.allpb = allpb;

    if ~perm %- keep more param when not perm
        ALL.TrialType = single(TrialType);
        ALL.trial2take = restr;
        ALL.T = T;
    end
    
    param.modelspec = modelspec;

    %- plot a quick recap figure, if needed
    
    if showfig
        
        figure;
        colors = cbrewer('div', 'PiYG', 64);
        colors = flipud(colors); % puts pink on top, green at the bottom

        subplot(3,4,1);
        imagesc(newf,[0 1]);axis xy
        colormap(colors);
        hold on
        plot(thr(:,1),1:length(probas),'-k','LineWidth',2);ylabel('Proba J1');xlabel('Proba J2');
        
        subplot(3,4,2)
        imagesc(conf);axis xy
        colormap(colors);
        title([mdl.Formula.LinearPredictor ' / BIC: ' num2str(mdl.ModelCriterion.BIC) ' / Disp= ' num2str(mdl.Dispersion) ' / Pref= ' num2str(pref) ' / R2= ' num2str(mdl.Rsquared.Adjusted)  ])
        
        subplot(3,4,[3 4]);
        bar(mdl.Coefficients.Estimate,'FaceColor',[.6 .6 .6])
        set(gca,'XTick',1:length(mdl.Coefficients.Estimate),'XTickLabels',mdl.CoefficientNames,'XTickLabelRotation',20	)
        pval =  round(mdl.Coefficients.pValue*10000)/10000;
        for i = 1 : length(mdl.Coefficients.tStat)
            if pval(i)<0.01
                text(i,mdl.Coefficients.Estimate(i),num2str(pval(i)),'Color','r','HorizontalAlignment','center')
            else
                text(i,mdl.Coefficients.Estimate(i),num2str(pval(i)),'HorizontalAlignment','center')
            end
        end
        xlim([0.5 length(pval)+.5])
        
        subplot(3,4,[5 6]);
        if p_trend<0.01 ; plot(mdl.Residuals.Raw,'r');
        else plot(mdl.Residuals.Raw,'k');
        end
        text(25,-0.75,['p-trend=' num2str(p_trend)])
        ylabel('Residuals');xlabel('Trials');
        
        subplot(3,4,[9 10]);
        line([0 length(pref_bins)],[0.5 0.5],'Color','k');box on;hold on
        if p_trend_bins<0.01 ; plot(pref_bins,'r')
        else plot(pref_bins,'k')
        end
        text(25,0.15,['p-trend=' num2str(p_trend_bins)])
        ylabel('Preference (J1 when > 0.5)');
        hold off
        ylim([0 1])
        
        subplot(3,4,[11]);
        if p_corr(1,2)<0.01
            plot(ft_bins,pref_bins,'.r')
        else
            plot(ft_bins,pref_bins,'.k')
        end
        xlabel('FT')
        title(['R=' num2str(r_corr(1,2)) ' / p=' num2str(p_corr(1,2))])
        ylim([0 1])
        
    end
    
else
    ALL = [];
    param = [];
end


