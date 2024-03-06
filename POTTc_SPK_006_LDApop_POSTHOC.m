%% LDA on pseudo-population - POST-PROCESSING - Figure 6
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.04

%% POST HOC - EACH MONKEY - FIGURE 6A-B-C

clear

nbNeur = 200;

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
%load([path2go 'res_LDApop_kfold_final_103sessions.mat'])
%load([path2go 'res_LDApop_kfold_final_186sessions.mat'])
load([path2go 'res_LDApop_kfold_final_289sessions.mat'])
measures = {'P_juice' 'P_proba' 'P_side' 'I_chosenjuice' 'I_chosenproba' 'I_unchosenproba' 'I_chosenside'}
measures = {'I_chosenjuice' 'I_chosenproba' 'I_unchosenproba' 'I_chosenside'}
evtCol = 'Greens';

%rew = load([path2go 'res_LDApop_rew_kfold_final_103sessions.mat'])
%rew = load([path2go 'res_LDApop_rew_kfold_final_186sessions.mat'])
rew = load([path2go 'res_LDApop_rew_kfold_final_289sessions.mat'])
allfields = fieldnames(rew.res_LDApop);
for i = 1 : length(allfields)
    res_LDApop.(allfields{i}) = rew.res_LDApop.(allfields{i});
end
measures = {'I_chosenjuice' 'I_chosenproba' 'I_chosenside' 'I_rew' 'I_chosenjuicerew'}
chance_l = [.5 .25 .5 .5 .5];
evtCol = {'Greens' 'Greens' 'Greens' 'Reds' 'Reds' };

colorsArea = cbrewer('qual', 'Paired', 12);
if length(area2test)==5
    colorsArea = colorsArea([2 4 5 6 12],:);
else
    colorsArea = colorsArea([1:7 12 8:11],:);
end

meanPerf={};
meanPerf_perm={};
nbNeur_sig =[];
for m = 1 : length(measures)
    x=0;
    for ar = 1: length(area2test)
        eval(['temp = res_LDApop.' measures{m} '.area_' area2test{ar} ';'])
        for n = 1 : length(temp)
            x = x+1;
            pval = sum(temp(n).perf_perm > nanmean(temp(n).perf)) / length(temp(n).perf_perm);
            sem = (nanstd(temp(n).perf)/sqrt(length(temp(n).perf)));
            meanPerf{m}(x,:) = [nanmean(temp(n).perf) temp(n).nbUnit ar nanmean(temp(n).perf)-sem nanmean(temp(n).perf)+sem pval];
            sem_perm = (nanstd(temp(n).perf_perm)/sqrt(length(temp(n).perf_perm)));
            meanPerf_perm{m}(x,:) = [nanmean(temp(n).perf_perm) temp(n).nbUnit ar nanmean(temp(n).perf_perm)-sem_perm nanmean(temp(n).perf_perm)+sem_perm];

        end
    end

    for ar = 1: length(area2test)
        idx = find(meanPerf{m}(:,3)==ar & meanPerf{m}(:,end)<0.05,1,'first');
        if ~isempty(idx)
            nbNeur_sig(m,ar) = meanPerf{m}(idx,2);
        else
            nbNeur_sig(m,ar) = NaN;
        end
        idx2 = find(meanPerf{m}(:,3)==ar ,1,'last');
        nbNeur_test(m,ar) = meanPerf{m}(idx2,2);
    end
end

%- Figure 6A,B,C
all_decoding_perf = 0;
for  m = 1  : length(measures)
    % fig = figure('Position',[1051 51 375 918]);
    fig = figure('Position',[1125 51 301 918]);
    subplot(3,1,[1 2]);
    for ar = 1: length(area2test)
        dumm = meanPerf{m}(meanPerf{m}(:,3)==ar,[1 2 4 5]);
        [param_a(m,ar),param_b(m,ar),Yh] = fitDecodPerf(dumm,chance_l(m));
        plot(dumm(:,2),dumm(:,1),'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15,'MarkerFaceColor',colorsArea(ar,:),'MarkerEdgeColor',colorsArea(ar,:));hold on
        plot(dumm(:,2),Yh+chance_l(m),'--','Color',colorsArea(ar,:),'LineWidth',1,'MarkerSize',15);hold on
        for p = 1:size(dumm,1)
            line([dumm(p,2) dumm(p,2)],[dumm(p,3) dumm(p,4)],'Color',colorsArea(ar,:))
        end
        keep4legend(ar,1) = dumm(1,1);

        dumm = meanPerf_perm{m}(meanPerf_perm{m}(:,3)==ar,[1 2 4 5]);
        plot(dumm(:,2),dumm(:,1),'.-','Color',colorsArea(ar,:),'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',colorsArea(ar,:),'MarkerEdgeColor',colorsArea(ar,:));hold on
        for p = 1:size(dumm,1)
            line([dumm(p,2) dumm(p,2)],[dumm(p,3) dumm(p,4)],'Color',colorsArea(ar,:))
        end

        % text(825,1-(ar/60),area2test{ar},'Color',colorsArea(ar,:),'FontWeight','bold','FontSize',16)
        %  ylim([0.4 1])
        xlim([-150 1100])

    end
    y_ax = min(keep4legend):(max(keep4legend)-min(keep4legend))/(length(area2test)-1) :max(keep4legend);
    [aa,bb]=sortrows(keep4legend);
    for ar = 1 : length(area2test)
        text(0,y_ax(ar),area2test{bb(ar)},'Color',colorsArea(bb(ar),:),'FontWeight','bold','FontSize',14,'HorizontalAlignment','right')
    end
    title(['Decoding perf - ' measures{m}])
    xlabel('Nb neurons included')
    ylabel('Average decoding perf')
    set(gca,'FontSize',16)

    decoding_perf = meanPerf{m}(meanPerf{m}(:,2)==nbNeur,[1 3]);

    area2plot = area2test(decoding_perf(:,2));
    for ar =1 : length(area2plot)
        if strcmp(area2plot{ar}(1),'1')
            area2plot{ar} = ['a' area2plot{ar}];
        end
    end

    if round(mean(meanPerf_perm{m}(:,1)),2)==0.5
        scales = [.5 .9];
    else
        scales = [.25 .7];
    end

    fig2 = subplot(3,1,3);
    brain_map(decoding_perf,area2plot,scales,evtCol{m},fig2)
    colorbar off
    colorbar('Location','South','Position',[0.337 0.121 0.459 0.009]);

    temp = meanPerf{m}(meanPerf{m}(:,2)==nbNeur,[1 3]);
    temp_perm = meanPerf_perm{m}(meanPerf_perm{m}(:,2)==nbNeur,[1 3]);

    all_decoding_perf = all_decoding_perf + (temp(:,1)-temp_perm(:,1))./(1-temp_perm(:,1));
end

% brain_map(all_decoding_perf/length(measures),area2plot)
% figure;
% [a,b]=sortrows(all_decoding_perf,"ascend")
% for ar = 1 : length(all_decoding_perf);
%     barh(ar,all_decoding_perf(b(ar))/length(measures),FaceColor=colorsArea(b(ar),:));hold on
% end
% set(gca,"YTick",1:length(area2test),'YTickLabel',area2test(b),'FontSize',16)
% xlabel('Information Content')

%% POST HOC - EACH MONKEY - Figure 6D

clear

nbNeur = 100;

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
stim_1 = load([path2go 'res_LDApop_kfold_final_103sessions.mat'])
stim_2 = load([path2go 'res_LDApop_kfold_final_186sessions.mat'])
stim_all = load([path2go 'res_LDApop_kfold_final_289sessions.mat'])
res_LDApop(1) = stim_1.res_LDApop;
res_LDApop(2) = stim_2.res_LDApop;
res_LDApop(3) = stim_all.res_LDApop;
measures = {'P_juice' 'P_proba' 'P_side' 'I_chosenjuice' 'I_chosenproba' 'I_unchosenproba' 'I_chosenside'}
measures = {'I_chosenjuice' 'I_chosenproba' 'I_unchosenproba' 'I_chosenside'}
evtCol = 'Greens';
area2test = stim_1.area2test;

rew_1 = load([path2go 'res_LDApop_rew_kfold_final_103sessions.mat'])
rew_2 = load([path2go 'res_LDApop_rew_kfold_final_186sessions.mat'])
rew_all = load([path2go 'res_LDApop_rew_kfold_final_289sessions.mat'])
allfields = fieldnames(rew_1.res_LDApop);
for i = 1 : length(allfields)
    res_LDApop(1).(allfields{i}) = rew_1.res_LDApop.(allfields{i});
    res_LDApop(2).(allfields{i}) = rew_2.res_LDApop.(allfields{i});
    res_LDApop(3).(allfields{i}) = rew_all.res_LDApop.(allfields{i});
end
measures = {'I_chosenjuice' 'I_chosenproba' 'I_chosenside' 'I_rew' 'I_chosenjuicerew'}
evtCol = {'Greens' 'Greens' 'Greens' 'Reds' 'Reds' };

colorsArea = cbrewer('qual', 'Paired', 12);
if length(area2test)==5
    colorsArea = colorsArea([2 4 5 6 12],:);
else
    colorsArea = colorsArea([1:7 12 8:11],:);
end

meanPerf={};
allPerf=cell(length(measures),1);
meanPerf_perm={};
nbNeur_sig =[];
for m = 1 : length(measures)
    x=0;
    for mk = 1 : 3
        for ar = 1: length(stim_1.area2test)
            eval(['temp = res_LDApop(mk).' measures{m} '.area_' stim_1.area2test{ar} ';'])
            for n = 1 : length(temp)
                x = x+1;
                pval = sum(temp(n).perf_perm > nanmean(temp(n).perf)) / length(temp(n).perf_perm);
                sem = (nanstd(temp(n).perf)/sqrt(length(temp(n).perf)));
                meanPerf{m}(x,:) = [nanmean(temp(n).perf) temp(n).nbUnit ar nanmean(temp(n).perf)-sem nanmean(temp(n).perf)+sem pval mk];
                sem_perm = (nanstd(temp(n).perf_perm)/sqrt(length(temp(n).perf_perm)));
                meanPerf_perm{m}(x,:) = [nanmean(temp(n).perf_perm) temp(n).nbUnit ar nanmean(temp(n).perf_perm)-sem_perm nanmean(temp(n).perf_perm)+sem_perm mk];
    
                nbRep = length(temp(n).perf');
                allPerf{m} = [allPerf{m} ; temp(n).perf' repmat(temp(n).nbUnit,nbRep,1) repmat(ar,nbRep,1) repmat(mk,nbRep,1)];
           end
        end
    end
end

tab_all={};
for m = 1 : length(measures)
    data_sub = meanPerf{m}(meanPerf{m}(:,end)~=3,[1 2 3 end]);
    mks = {'M' , 'X'}  ;
    %- significance
    modeldata_pop = table(data_sub(:,1),data_sub(:,2),stim_1.area2test(data_sub(:,3))',mks(data_sub(:,4))','VariableNames',{'perf' 'nbunit' 'area' 'mk'});
    models_form = {'perf ~ 1 + area  + (1|mk) + (1|nbunit)' ; 'perf ~ 1 + area  + (1|mk)'};
    % [lme,model_final] = model_comparison(modeldata_pop,models_form,false);
    lme = fitglme(modeldata_pop,models_form{1}); model_final = models_form{1};
    
    [pval,wald,thr_corr,pval_adj] = area_posthoc(lme,stim_1.area2test,'y');
    disp(['%%%%%%%%%% Decoding perf with model = ' model_final ' %%%%%%%%%%'])
    disp(anova(lme));disp(pval_adj);disp(wald);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    [tab, arlist] = build_table(pval_adj,wald,area2test);
    tab_all = [tab_all , tab];
    pause
end

all_decoding_perf = 0;
fig = figure;
for  m = 1  : length(measures)
    % fig = figure('Position',[1051 51 375 918]);
    subplot(1,5,m);
    decoding_perf = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==3,[1 4 5 3]);
    decoding_perf_M = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==1,[1 4 5 3]);
    decoding_perf_X = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==2,[1 4 5 3]);
    area2plot = stim_1.area2test(decoding_perf(:,end));
    for ar = 1 : length(area2plot)
        if strcmp(area2plot{ar}(1),'1')
            area2plot{ar} = ['a' area2plot{ar}];
        end
    end

    for ar = 1 : length(area2plot)
        bar(decoding_perf(ar,end),decoding_perf(ar,1),'FaceColor',colorsArea(ar,:));hold on
        line([ar ar], [decoding_perf(ar,2) decoding_perf(ar,3) ],'Color','k');hold on
    end
    plot(decoding_perf_M(:,end),decoding_perf_M(:,1),'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot(decoding_perf_X(:,end),decoding_perf_X(:,1),'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    ylim([0 1])
    xlim([0 length(area2plot)+1])
    set(gca,'Xtick',1:length(area2plot),'XtickLabel',area2plot,'XtickLabelRotation',30,'FontSize',16)
    title(measures{m})
end
