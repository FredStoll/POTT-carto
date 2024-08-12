%% POST HOC FOR ANOVA + LDA = FIGURES 3-4-5-8 

%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.03

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'

period = 'STIM' ; % STIM

%- load the anova results
load([path2go 'res_ANOVA_full_nested.mat'])

if strcmp(period,'REW')
    params = [1 2];
    mat = 'rew_diff';
    measures = {'I_rew' 'I_chosenjuice_rew'};
    LDA = load([path2go 'res_LDA_rew_4neurons.mat']);

    %- periods for computing average number of sig neurons during time bins : name / time alignment / times
    % periods = {'REF'  , 'Stim_onset' , [-600 -100]  ;
    periods = {'REF'  , 'FixFP_onset' , [100 700]  ;
               'REW'  ,  'Rew_onset' ,  [100 700] };
    %- second column refers to times_evts matrix, like: times_evts = {'FixFP_onset' 'Stim_onset' 'Resp_onset' 'FixResp_onset' 'FB_onset' 'Rew_onset' 'FB_offset'};

else
     params = [2 1 ];
     measures = {'I_chosenjuice' 'I_chosenproba'};
%    params = [3 1];
%    measures = {'I_chosenside' 'I_chosenproba'};
    mat = 'ins_diff';
    LDA = load([path2go 'res_LDA_stim_4neurons.mat']);

    %- periods for computing average number of sig neurons during time bins : name / time alignment / times
    % periods = {'REF'  , 'Stim_onset' , [-600 -100]  ;
    periods = {'REF'  , 'FixFP_onset' , [100 700]  ;
                'STIM' , 'Stim_onset' ,  [100 700] };
    %- second column refers to times_evts matrix, like: times_evts = {'FixFP_onset' 'Stim_onset' 'Resp_onset' 'FixResp_onset' 'FB_onset' 'Rew_onset' 'FB_offset'};
end
periods_all = {'REF'  , 'FixFP_onset' , [100 700]  ;
               'STIM' , 'Stim_onset' ,  [100 700] ;
               'REW'  ,  'Rew_onset' ,  [100 700] };

thr=0.01;nsig=3; %- threshold for sig

%- for plotting over time
step = 20; % step in time for plotting purposes only (write down value every X time bin)
timesel = [2 4 5 6];
timesel_sub = [-.4 .98 ; 0 .98 ; -.2 .38 ; -.2 .8];
areas = utils_POTT_areas;
area2test = {'a12r' 'a12m' 'a12o' 'a12l' 'a11ml' 'a13l' 'a13m' 'LAI' };

%- COLOR ASSIGNMENT
colorsArea = cbrewer('qual', 'Paired', 12);
if length(area2test)==5
    colorsArea = colorsArea([2 4 5 6 12],:);
else
    colorsArea = colorsArea([1:7 12 8:11],:);
end
colorsArea_sub = colorsArea;

x = 0;
for ar = 1 : length(area2test)
    eval(['take = LDA.res_LDA.I_chosenjuice.' area2test{ar} ';'])
    for i = 1 : length(take)
        x=x+1;
        nUnits(x) = length(take(i).takeit);
    end
end 
prctile(nUnits,[50 0 100])

all_diff = [];
all_diff_both = [];
all_diff_omega = [];
all_diff_betas = [];
all_diff_betas_both = [];
all_sess = [];
all_units = [];
all_mk = [];
all_neuron_id = [];
all_neuron_info = [];
for s = 1 : length(res_anova)
    if ~isempty(res_anova(s).ins_diff) & ~isempty(res_anova(s).rew_diff)
        all_diff_omega = [all_diff_omega ; res_anova(s).(mat).Omega2];
        all_diff_betas = [all_diff_betas ; res_anova(s).(mat).Betas];
        all_diff = [all_diff ; res_anova(s).(mat).Ap];
        all_diff_both = [all_diff_both ; cat(2,res_anova(s).ins_diff.Ap,res_anova(s).rew_diff.Ap)];
        all_diff_betas_both = [all_diff_betas_both ; cat(2,res_anova(s).ins_diff.Betas,res_anova(s).rew_diff.Betas)];

        all_units = [all_units ; res_anova(s).neurons_area];
        all_mk = [all_mk ; repmat({res_anova(s).session(1)},size(res_anova(s).neurons_area,1),1)];
        all_sess = [all_sess ; repmat(s,size(res_anova(s).neurons_area,1),1)];

        for n = 1 : length(res_anova(s).neurons_area)
            %- extract the name of the file where the timestamps are!
            if res_anova(s).neurons_info(n,1)<10; addzeros = '00';
            elseif res_anova(s).neurons_info(n,1)<100; addzeros = '0';
            else addzeros = '';
            end
            if res_anova(s).neurons_info(n,2)<10; addzeros2 = '0';
            else addzeros2 = '';
            end

            all_neuron_id = [all_neuron_id ; [res_anova(s).session '_Ch' addzeros num2str(res_anova(s).neurons_info(n,1)) ...
                '_Clus' addzeros2 num2str(res_anova(s).neurons_info(n,2)) ]];
            all_neuron_info = [all_neuron_info; res_anova(s).neurons_info(n,1) res_anova(s).neurons_info(n,3)];

        end

    end
end

for p = 1 : size(all_diff,2)
    dumm = squeeze(all_diff(:,p,:));
    pval_sig = zeros(size(dumm));
    for k = 1 : length(dumm(:,1)) %- for each neuron, check if below thr
        [~,idxs] = findenough(dumm(k,:),thr,nsig,'<=');
        pval_sig(k,idxs)=1;clear idxs
    end
    all_diff_sig(:,p,:) = pval_sig;
end

for p = 1 : size(all_diff_both,2)
    dumm = squeeze(all_diff_both(:,p,:));
    pval_sig = zeros(size(dumm));
    for k = 1 : length(dumm(:,1)) %- for each neuron, check if below thr
        [~,idxs] = findenough(dumm(k,:),thr,nsig,'<=');
        pval_sig(k,idxs)=1;clear idxs
    end
    all_diff_both_sig(:,p,:) = pval_sig;
end


%- remove neurons where models fails
keep = (sum(isnan(squeeze(all_diff_both(:,1,:)))')==0)' & (sum(isnan(squeeze(all_diff_both(:,end,:)))')==0)';
mk1 = ismember(all_mk,'M')

all_units_area = cell(size(all_units));
for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,areas.' area2test{ar} ');'])
    all_units_area(takeit) = area2test(ar);
end

%- remove unlabeled neurons
keep = keep & ~cellfun(@isempty,all_units_area);

%- number of units // sessions (WHEN ENOUGH TRIALS FOR FIG2!!!)
 for ar = 1 : length(area2test)
     eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
     nbunit_tot(ar,:) = [sum(takeit & mk1) sum(takeit & ~mk1)];
  %   nbsess_tot(ar,:) = [length(unique(all_sess(takeit & mk1))) length(unique(all_sess(takeit & ~mk1)))];
 end

%% reactivation - Fig 8B

%periods_reac = periods_all(2:3,:);
periods_reac = {'STIM' , 'Stim_onset' ,  [100 700] ;
                  'HOLD'  ,   'FixResp_onset' ,  [100 600]  ;
                  'FB'  ,   'FB_onset' ,  [0 500]  ;
           %     'REW'  ,  'Rew_onset' ,  [100 700] };
                 'REW'  ,  'Rew_onset' ,  [100 600] };
periods2compare = [1 2 ; 1 3 ; 1 4];
%periods2compare = [3 4];
params_reac = [2 1 3];
xx=0;
prop_sig_all=[];
clear PVAL_all RVAL_all allstats
figure
for p = 1 : length(params_reac)
    for pe = 1 : length(periods2compare(:,1))
        time_chunk = bins_considered == find(ismember(times_evts,periods_reac{periods2compare(pe,1),2}));
        time_considered = find(time_chunk & time>=periods_reac{periods2compare(pe,1),3}(1) & time<=periods_reac{periods2compare(pe,1),3}(2));

        time_chunk2 = bins_considered == find(ismember(times_evts,periods_reac{periods2compare(pe,2),2}));
        time_considered2 = find(time_chunk2 & time>=periods_reac{periods2compare(pe,2),3}(1) & time<=periods_reac{periods2compare(pe,2),3}(2));

        nb_sig_reac = [];
        for ar = 1 : length(area2test)
            eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
            dumm = squeeze(all_diff_sig(:,params_reac(p),:))';
            sig_reac = [(sum(dumm(time_considered,takeit))~=0)' (sum(dumm(time_considered2,takeit))~=0)'];
            sig_reac_M = [(sum(dumm(time_considered,takeit & mk1))~=0)' (sum(dumm(time_considered2,takeit & mk1))~=0)'];
            sig_reac_X = [(sum(dumm(time_considered,takeit & ~mk1))~=0)' (sum(dumm(time_considered2,takeit & ~mk1))~=0)'];

            %         nb_sig_reac(ar,:) = [sum(sig_reac(:,1)==1 & sig_reac(:,2)==0) ...
            %                              sum(sig_reac(:,1)==1 & sig_reac(:,2)==1) ...
            %                              sum(sig_reac(:,1)==0 & sig_reac(:,2)==1) ...
            %                              sum(sig_reac(:,1)==0 & sig_reac(:,2)==0) ];
            nb_sig_reac(ar,:) = [sum(sig_reac(:,1)==1 & sig_reac(:,2)==0) ...
                sum(sig_reac(:,1)==1 & sig_reac(:,2)==1) ];
            %                      sum(sig_reac(:,1)==0 & sig_reac(:,2)==0) ];
            %   nb_sig_reac_M(ar,:) = [sum(sig_reac_M(:,1)) sum(sum(sig_reac_M')==2) sum(sig_reac_M(:,2)) size(sig_reac_M,1)];
            %   nb_sig_reac_X(ar,:) = [sum(sig_reac_X(:,1)) sum(sum(sig_reac_X')==2) sum(sig_reac_X(:,2)) size(sig_reac_X,1)];

            %- extract average beta param during the different windows
            betas_temp = squeeze(all_diff_betas(takeit,params_reac(p),:));
            names_temp = all_neuron_id(takeit,:);
            beta_t1 = nanmean(betas_temp(sig_reac(:,1)==1 & sig_reac(:,2)==1,time_considered)');
            beta_t2 = nanmean(betas_temp(sig_reac(:,1)==1 & sig_reac(:,2)==1,time_considered2)');
            signeurons_id = names_temp(sig_reac(:,1)==1 & sig_reac(:,2)==1,:);
            [rval,pval] = corrcoef(beta_t1,beta_t2);
            PVAL_all{p}(ar,pe) = pval(2,1);
            RVAL_all{p}(ar,pe) = rval(2,1);
            if ar == 4 & p == 2
                xx = xx + 1;
                example_betas = [beta_t1;beta_t2];
                subplot(1,3,xx);
                plot(example_betas(1,:),example_betas(2,:),'.r');hold on

                %- find my example neuron!
                takethat = ismember(signeurons_id,{'X032621a_Ch077_Clus01'});
                plot(example_betas(1,takethat),example_betas(2,takethat),'.k');hold on

                cf = fit(example_betas(1,:)',example_betas(2,:)','poly1'); % fit
                p_handle = plot(cf,'k','predfunc',.95); % Plot fit
                set(p_handle,'Color','k','LineWidth',2);
                disp(RVAL_all{p}(ar,pe))
                disp(PVAL_all{p}(ar,pe))
                xlim([-.25 .25])
                ylim([-.25 .25])
                title([area2test{ar} ' - ' res_anova(1).ins_diff.predic{params_reac(p)}])
            end
        end

        prop_sig= nb_sig_reac./repmat(sum(nb_sig_reac')',1,size(nb_sig_reac,2));
        prop_sig_all(p,pe,:) = prop_sig(:,2)';

    end
end

%- for the stats
for p = 1 : length(params_reac)
    for pe = 1 : length(periods2compare(:,1))
        time_chunk = bins_considered == find(ismember(times_evts,periods_reac{periods2compare(pe,1),2}));
        time_considered = find(time_chunk & time>=periods_reac{periods2compare(pe,1),3}(1) & time<=periods_reac{periods2compare(pe,1),3}(2));

        time_chunk2 = bins_considered == find(ismember(times_evts,periods_reac{periods2compare(pe,2),2}));
        time_considered2 = find(time_chunk2 & time>=periods_reac{periods2compare(pe,2),3}(1) & time<=periods_reac{periods2compare(pe,2),3}(2));

        dumm = squeeze(all_diff_sig(:,params_reac(p),:))';
        sig_reac_all = [(sum(dumm(time_considered,:))~=0)' (sum(dumm(time_considered2,:))~=0)'];
        sig = NaN(length(keep),1);
        sig(sig_reac_all(:,1)==1 & sig_reac_all(:,2)==0,1)=false;
        sig(sig_reac_all(:,1)==1 & sig_reac_all(:,2)==1,1)=true;

        modeldata_su = table(sig(keep),all_units_area(keep),all_mk(keep),categorical(all_sess(keep)), 'VariableNames',{'sig' 'area' 'mk' 'sess'});
        modeldata_su(isnan(modeldata_su.sig),:)=[];
        modeldata_su.sig = (modeldata_su.sig)==1;

        models_form = {'sig ~ 1 + area  + (1|mk) + (1|sess)'};
        lme = fitglme(modeldata_su,models_form{1},'Distribution','Binomial');
        [~,wald_percNeurons,~,pval_adj_percNeurons]  = area_posthoc(lme,area2test,'y');
         title(['%%%%%%%%%% Percent of ' res_anova(1).(mat).predic{params_reac(p)} ' neurons during ' num2str(periods2compare(pe,1)) ' vs ' num2str(periods2compare(pe,2))  ' %%%%%%%%%%'])
       disp(['%%%%%%%%%% Percent of ' res_anova(1).(mat).predic{params_reac(p)} ' neurons during ' num2str(periods2compare(pe,1)) ' vs ' num2str(periods2compare(pe,2))  ' %%%%%%%%%%'])
        disp(anova(lme));disp(pval_adj_percNeurons);disp(wald_percNeurons);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        allstats{p,pe} = pval_adj_percNeurons;
    end
end

figure;
colorsPeriods = [158 233 255 ; 253 176 176 ; 255 0 0]/255;
for p = 1 : length(params_reac)
    subplot(1,3,p);
    data2show = squeeze(prop_sig_all(p,:,:))';
    bb = bar(squeeze(prop_sig_all(p,:,:))');
    for b = 1 : length(bb)
        set(bb(b),'FaceColor',colorsPeriods(b,:),'DisplayName',periods_reac{periods2compare(b,2),1})
    end
    legend
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'FontSize',16)
    title(res_anova(1).ins_diff.predic{params_reac(p)}  )
    ylim([0 1])
    xlab = [-0.25 0 0.25]-0.1
    old=0.01*ones(length(area2test),3);
    for pe = 1:size(allstats,2)
        temp_stats =  allstats{p,pe};
        temp_stats = [temp_stats zeros(8,1)] + [temp_stats zeros(8,1)]';
        temp_stats(temp_stats==0)=NaN;
        for a1 = 1 : size(temp_stats,1)
            for a2 = 1 : size(temp_stats,2)
                if temp_stats(a1,a2)<0.05 & data2show(a1,pe)>data2show(a2,pe)
                    text(a1+xlab(pe),data2show(a1,pe)+old(a1,pe),'*','Color',colorsArea(a2,:),'FontSize',30)
                    old(a1,pe) = old(a1,pe)+0.0175;
                end
            end
        end
    end
end

figure;
colors = cbrewer('seq', 'Greys', 100); 
rand_x = [-.3:.6/7:.3];
for p = 1 : length(params_reac)
    subplot(3,1,p);
    data2show = squeeze(prop_sig_all(p,:,:));
    imagesc(data2show,[0 1]);colormap(colors)
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'FontSize',16)
    set(gca,'Ytick',1:3,'YtickLabel',{'t(Resp)' 't(FB)' 't(Rew)'},'FontSize',16)
    title(res_anova(1).ins_diff.predic{params_reac(p)}  )
    for pe = 1:size(allstats,2)
        temp_stats =  allstats{p,pe};
        temp_stats = [temp_stats zeros(8,1)] + [temp_stats zeros(8,1)]';
        temp_stats(temp_stats==0)=NaN;
        for a1 = 1 : size(temp_stats,1)
            for a2 = 1 : size(temp_stats,2)
                if temp_stats(a1,a2)<0.05 & data2show(pe,a1)>data2show(pe,a2)
                    text(a1+rand_x(a2),pe,'*','Color',colorsArea(a2,:),'FontSize',24)
                end
            end
        end
    end
end

figure;
colorsPeriods = [158 233 255 ; 253 176 176 ; 255 0 0]/255;
for p = 1 : length(params_reac)
    subplot(1,3,p);
    bb = bar(RVAL_all{p}); hold on
    for b = 1 : length(bb)
        set(bb(b),'FaceColor',colorsPeriods(b,:),'DisplayName',periods_reac{periods2compare(b,2),1})
    end
    legend
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'FontSize',16)
    title(res_anova(1).ins_diff.predic{params_reac(p)}  )
    ylim([-.25 1])

    [h, crit_p, adj_p]=fdr_bh(PVAL_all{p},.05);
    xax = [-0.25 0 0.25];
    for ar = 1 : length(area2test)
        for t = 1 : size(periods2compare,1)
            text(ar+xax(t),-0.2,num2str(adj_p(ar,t)),Rotation=90)
        end
    end

end

figure;
colors = cbrewer('seq', 'Greys', 100);
rand_x = [-.3:.6/7:.3];
for p = 1 : length(params_reac)
    subplot(3,1,p);
    imagesc(RVAL_all{p}',[-.2 1]);colormap(colors)
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'FontSize',16)
    set(gca,'Ytick',1:3,'YtickLabel',{'t(Resp)' 't(FB)' 't(Rew)'},'FontSize',16)
    title(res_anova(1).ins_diff.predic{params_reac(p)}  )
    [h, crit_p, adj_p]=fdr_bh(PVAL_all{p},.05);

    for pe = 1:size(adj_p,2)

        for a1 = 1 : size(adj_p,1)
            if adj_p(a1,pe)<0.05
                text(a1,pe,'*','Color','k','FontSize',24)
            end
        end
    end
end

%% Fig 5D - venn for chosen flavor neurons at reward

params = [1 2];
mat2take = 'rew_diff';
measures2take = {'I_rew' 'I_chosenjuice_rew'};


%periods_reac = periods_all(2:3,:);
periods_reac = {'STIM' , 'Stim_onset' ,  [100 700] ;
                'HOLD'  ,   'FixResp_onset' ,  [100 600]  ;
                'FB'  ,   'FB_onset' ,  [0 500]  ;
                %     'REW'  ,  'Rew_onset' ,  [100 700] };
                'REW'  ,  'Rew_onset' ,  [100 600] };
periods2compare = [4 1];
params_reac = [2];
x = 0;
prop_sig_all=[];
for p = 1 : length(params_reac)
    for pe = 1 : length(periods2compare(:,1))
        x = x + 1;
        time_chunk = bins_considered == find(ismember(times_evts,periods_reac{periods2compare(pe,1),2}));
        time_considered = find(time_chunk & time>=periods_reac{periods2compare(pe,1),3}(1) & time<=periods_reac{periods2compare(pe,1),3}(2));

        time_chunk2 = bins_considered == find(ismember(times_evts,periods_reac{periods2compare(pe,2),2}));
        time_considered2 = find(time_chunk2 & time>=periods_reac{periods2compare(pe,2),3}(1) & time<=periods_reac{periods2compare(pe,2),3}(2));

        nb_sig_reac = [];nb_sig_reac_M=[];nb_sig_reac_X=[]
        for ar = 1 : length(area2test)
            eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
            dumm = squeeze(all_diff_sig(:,params_reac(p),:))';
            sig_reac = [(sum(dumm(time_considered,takeit))~=0)' (sum(dumm(time_considered2,takeit))~=0)'];
            sig_reac_M = [(sum(dumm(time_considered,takeit & mk1))~=0)' (sum(dumm(time_considered2,takeit & mk1))~=0)'];
            sig_reac_X = [(sum(dumm(time_considered,takeit & ~mk1))~=0)' (sum(dumm(time_considered2,takeit & ~mk1))~=0)'];

            nb_sig_reac(ar,:) = [sum(sig_reac(:,1)==1 & sig_reac(:,2)==0) ...
                                sum(sig_reac(:,1)==1 & sig_reac(:,2)==1) ...
                                sum(sig_reac(:,1)==0 & sig_reac(:,2)==1) ...
                                sum(sig_reac(:,1)==0 & sig_reac(:,2)==0) ];
            nb_sig_reac_M(ar,:) = [sum(sig_reac_M(:,1)==1 & sig_reac_M(:,2)==0) ...
                                sum(sig_reac_M(:,1)==1 & sig_reac_M(:,2)==1) ...
                                sum(sig_reac_M(:,1)==0 & sig_reac_M(:,2)==1) ...
                                sum(sig_reac_M(:,1)==0 & sig_reac_M(:,2)==0) ];
            nb_sig_reac_X(ar,:) = [sum(sig_reac_X(:,1)==1 & sig_reac_X(:,2)==0) ...
                                sum(sig_reac_X(:,1)==1 & sig_reac_X(:,2)==1) ...
                                sum(sig_reac_X(:,1)==0 & sig_reac_X(:,2)==1) ...
                                sum(sig_reac_X(:,1)==0 & sig_reac_X(:,2)==0) ];
        end

        prop_sig= nb_sig_reac./repmat(sum(nb_sig_reac')',1,size(nb_sig_reac,2));

        figure;
        for ar = 1 : length(area2test)
            subplot(1,length(area2test),ar)
            [H,S] = venn((nb_sig_reac(ar,[1 3])+nb_sig_reac(ar,2))/sum(nb_sig_reac(ar,:)), nb_sig_reac(ar,2)/sum(nb_sig_reac(ar,:)));
            for i = 1:3
                text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), [num2str(nb_sig_reac(ar,i))])
            end
            xlim([-.4 .7]);ylim([-.4 .4])
        end

        prop_sig_all(p,pe,:) = prop_sig(:,2)';
    end
end

%- neurons encoding 1 param vs 2+
nN = nb_sig_reac(:,1:3)./sum(nb_sig_reac(:,1:3),2);
nN_M = nb_sig_reac_M(:,1:3)./sum(nb_sig_reac_M(:,1:3),2);
nN_X = nb_sig_reac_X(:,1:3)./sum(nb_sig_reac_X(:,1:3),2);

figure;
bb=barh(nN,'stacked');axis ij;hold on
bb(1).FaceColor=[.6 .6 .6];
bb(2).FaceColor=[92 161 92]/255;
bb(3).FaceColor=[130 200 130]/255;
plot(nN_M(:,1) ,0.2+(1:length(area2test)),'o','MarkerSize',10,'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor','k')
plot(nN_M(:,1)+nN_M(:,2) ,0.2+(1:length(area2test)),'o','MarkerSize',10,'MarkerFaceColor',[130 200 130]/255,'MarkerEdgeColor','k')
plot(nN_X(:,1)  ,-0.2+(1:length(area2test)),'v','MarkerSize',10,'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor','k')
plot(nN_X(:,1)+nN_X(:,2) ,-0.2+(1:length(area2test)),'v','MarkerSize',10,'MarkerFaceColor',[130 200 130]/255,'MarkerEdgeColor','k')
set(gca,'Ytick',1:length(area2test),'YtickLabel',area2test,'FontSize',16)
xlabel('Proportion of neurons')

%% mixed selectivity across a few electrodes covering 12o/12l in mimic

% elecloc = load('C:\Users\Fred\Dropbox\Rudebeck Lab\SCRIPTS\POTT_RECORDINGS\Mimic_ElecLocations_Drive2_postHisto.mat')
% 
% params2take = [2 1 3];
% mat2take = 'ins_diff';
% 
% %periods_reac = periods_all(2:3,:);
% periods_mixed = {'STIM' , 'Stim_onset' ,  [100 700]  };
% 
% time_chunk = bins_considered == find(ismember(times_evts,periods_mixed{1,2}));
% time_considered = find(time_chunk & time>=periods_mixed{1,3}(1) & time<=periods_mixed{1,3}(2));
% elec2test = [55 56 66 67 68 69 70 77 78];
% figure;
% nb_sig_reac = [];
% for ar = 1 : length(elec2test)
%     sig_reac = [];
% 
%     takeit = all_neuron_info(:,1)==elec2test(ar) & keep & ~mk1 ;
%     depth_info = all_neuron_info(takeit,2);
% 
%     for p = 1 : length(params2take)
%         dumm = squeeze(all_diff_sig(:,params2take(p),:))';
%         sig_reac(:,p) = (sum(dumm(time_considered,takeit))~=0)';
%     end
%    
%     sg_param = (sig_reac(:,1)==1 & sig_reac(:,2)==0 & sig_reac(:,3)==0) | (sig_reac(:,1)==0 & sig_reac(:,2)==1 & sig_reac(:,3)==0)  | (sig_reac(:,1)==0 & sig_reac(:,2)==0 & sig_reac(:,3)==1) ;
%     mult_param = (sig_reac(:,1)==1 & sig_reac(:,2)==1 & sig_reac(:,3)==1) | (sig_reac(:,1)==1 & sig_reac(:,2)==1 & sig_reac(:,3)==0) | (sig_reac(:,1)==1 & sig_reac(:,2)==0 & sig_reac(:,3)==1)  | (sig_reac(:,1)==0 & sig_reac(:,2)==1 & sig_reac(:,3)==1);
% 
%     plot(elec2test(ar)-0.25+(rand(sum(sg_param),1)/4),depth_info(sg_param),'.k');hold on
%     plot(elec2test(ar)+0.25+(rand(sum(mult_param),1)/4),depth_info(mult_param),'.r');
% 
%     %- plot the grey matter
%     [lin,col]=find(elecloc.ELECSnb==elec2test(ar));
%     grey = -elecloc.ELECS{lin,col};
%     areanames = elecloc.AREA(elec2test(ar),:);
% 
%     for ii = 1 : size(grey,1)
%         if ismember(areanames{ii},{'12l'})
%             colorMe = colorsArea(4,:);
%         elseif ismember(areanames{ii},{'12o'})
%             colorMe = colorsArea(3,:);
%         else
%             colorMe='k'
%         end
%         line([elec2test(ar) elec2test(ar)],[grey(ii,1) grey(ii,2)],'Color',colorMe,'LineWidth',2)
%         line([elec2test(ar)-0.2 elec2test(ar)+.2],[grey(ii,1) grey(ii,1)],'Color',colorMe)
%         line([elec2test(ar)-0.2 elec2test(ar)+.2],[grey(ii,2) grey(ii,2)],'Color',colorMe)
%     end
% 
% end

%% Fig 5A-B - mixed selectivity at the Stimulus onset 

params2take = [2 1 3];
mat2take = 'ins_diff';

%periods_reac = periods_all(2:3,:);
periods_mixed = {'STIM' , 'Stim_onset' ,  [100 700]  };

time_chunk = bins_considered == find(ismember(times_evts,periods_mixed{1,2}));
time_considered = find(time_chunk & time>=periods_mixed{1,3}(1) & time<=periods_mixed{1,3}(2));

nb_sig_reac = [];nb_sig_reac_M=[];nb_sig_reac_X=[];
for ar = 1 : length(area2test)
    sig_reac = [];
    sig_reac_M = [];
    sig_reac_X = [];
    eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
    for p = 1 : length(params2take)
        dumm = squeeze(all_diff_sig(:,params2take(p),:))';
        sig_reac(:,p) = (sum(dumm(time_considered,takeit))~=0)';
        sig_reac_M(:,p) = (sum(dumm(time_considered,takeit & mk1))~=0)';
        sig_reac_X(:,p) = (sum(dumm(time_considered,takeit & ~mk1))~=0)';
    end

    % [z1 z2 z3 z12 z13 z23 z123]
    nb_sig_reac(ar,:) = [sum(sig_reac(:,1)==1 & sig_reac(:,2)==0 & sig_reac(:,3)==0) ...
        sum(sig_reac(:,1)==0 & sig_reac(:,2)==1 & sig_reac(:,3)==0) ...
        sum(sig_reac(:,1)==0 & sig_reac(:,2)==0 & sig_reac(:,3)==1) ...
        sum(sig_reac(:,1)==1 & sig_reac(:,2)==1 & sig_reac(:,3)==0) ...
        sum(sig_reac(:,1)==1 & sig_reac(:,2)==0 & sig_reac(:,3)==1) ...
        sum(sig_reac(:,1)==0 & sig_reac(:,2)==1 & sig_reac(:,3)==1) ...
        sum(sig_reac(:,1)==1 & sig_reac(:,2)==1 & sig_reac(:,3)==1) ...
        sum(sig_reac(:,1)==0 & sig_reac(:,2)==0 & sig_reac(:,3)==0) ];

    nb_sig_reac_M(ar,:) = [sum(sig_reac_M(:,1)==1 & sig_reac_M(:,2)==0 & sig_reac_M(:,3)==0) ...
        sum(sig_reac_M(:,1)==0 & sig_reac_M(:,2)==1 & sig_reac_M(:,3)==0) ...
        sum(sig_reac_M(:,1)==0 & sig_reac_M(:,2)==0 & sig_reac_M(:,3)==1) ...
        sum(sig_reac_M(:,1)==1 & sig_reac_M(:,2)==1 & sig_reac_M(:,3)==0) ...
        sum(sig_reac_M(:,1)==1 & sig_reac_M(:,2)==0 & sig_reac_M(:,3)==1) ...
        sum(sig_reac_M(:,1)==0 & sig_reac_M(:,2)==1 & sig_reac_M(:,3)==1) ...
        sum(sig_reac_M(:,1)==1 & sig_reac_M(:,2)==1 & sig_reac_M(:,3)==1) ...
        sum(sig_reac_M(:,1)==0 & sig_reac_M(:,2)==0 & sig_reac_M(:,3)==0) ];

    nb_sig_reac_X(ar,:) = [sum(sig_reac_X(:,1)==1 & sig_reac_X(:,2)==0 & sig_reac_X(:,3)==0) ...
        sum(sig_reac_X(:,1)==0 & sig_reac_X(:,2)==1 & sig_reac_X(:,3)==0) ...
        sum(sig_reac_X(:,1)==0 & sig_reac_X(:,2)==0 & sig_reac_X(:,3)==1) ...
        sum(sig_reac_X(:,1)==1 & sig_reac_X(:,2)==1 & sig_reac_X(:,3)==0) ...
        sum(sig_reac_X(:,1)==1 & sig_reac_X(:,2)==0 & sig_reac_X(:,3)==1) ...
        sum(sig_reac_X(:,1)==0 & sig_reac_X(:,2)==1 & sig_reac_X(:,3)==1) ...
        sum(sig_reac_X(:,1)==1 & sig_reac_X(:,2)==1 & sig_reac_X(:,3)==1) ...
        sum(sig_reac_X(:,1)==0 & sig_reac_X(:,2)==0 & sig_reac_X(:,3)==0) ];

end

figure;
for ar = 1 : length(area2test)
    subplot(1,length(area2test),ar)
    [H,S] = venn(nb_sig_reac(ar,1:end-1)/sum(nb_sig_reac(ar,:)));
    % [H,S] = venn((nb_sig_reac(ar,[1 3])+nb_sig_reac(ar,2))/sum(nb_sig_reac(ar,1:3)), nb_sig_reac(ar,2)/sum(nb_sig_reac(ar,1:3)));
    % [H,S] = venn(nb_sig_reac(ar,[1 3])+nb_sig_reac(ar,2), nb_sig_reac(ar,2));
    for i = 1:7
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), [num2str(nb_sig_reac(ar,i))],'HorizontalAlignment','center')
    end
    xlim([-.6 .8]);ylim([-.6 .8])
end


%- neurons encoding 1 param vs 2+
nN = [sum(nb_sig_reac(:,1:3),2) sum(nb_sig_reac(:,4:7),2) ];
nN_M = [sum(nb_sig_reac_M(:,1:3),2) sum(nb_sig_reac_M(:,4:7),2) ];
nN_X = [sum(nb_sig_reac_X(:,1:3),2) sum(nb_sig_reac_X(:,4:7),2) ];
clear pval
for ar = 1 : length(area2test)
     [tbl,chi2stat(ar),pval(ar)] = chi2_fms(nN(ar,1),nN(ar,1)+nN(ar,2),nN(ar,2),nN(ar,1)+nN(ar,2));
end
[h, crit_p, adj_p]=fdr_bh(pval',.05);

figure;
subplot(2,1,1)
bb=barh(nN./sum(nN,2),'stacked');axis ij;hold on
bb(1).FaceColor=[.7 .7 .7];
bb(2).FaceColor=[.3 .3 .3];
plot(nN_M(:,1)./sum(nN_M,2),1:length(area2test),'o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(nN_X(:,1)./sum(nN_X,2),1:length(area2test),'v','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k')
for ar = 1 : length(area2test)
    text(0.05,ar,['Chi2=' num2str(chi2stat(ar)) ', p=' num2str(adj_p(ar))])
end
set(gca,'Ytick',1:length(area2test),'YtickLabel',area2test,'FontSize',16)
xlabel('Proportion of neurons')

%- neurons encoding 1 param vs 3
nN = [sum(nb_sig_reac(:,1:3),2) nb_sig_reac(:,7) ];
nN_M = [sum(nb_sig_reac_M(:,1:3),2) nb_sig_reac_M(:,7) ];
nN_X = [sum(nb_sig_reac_X(:,1:3),2) nb_sig_reac_X(:,7) ];
for ar = 1 : length(area2test)
     [tbl,chi2stat(ar),pval(ar)] = chi2_fms(nN(ar,1),nN(ar,1)+nN(ar,2),nN(ar,2),nN(ar,1)+nN(ar,2));
end
[h, crit_p, adj_p]=fdr_bh(pval',.05);


subplot(2,1,2);
bb=barh(nN./sum(nN,2),'stacked');axis ij;hold on
bb(1).FaceColor=[.7 .7 .7];
bb(2).FaceColor=[.3 .3 .3];
plot(nN_M(:,1)./sum(nN_M,2),1:length(area2test),'o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(nN_X(:,1)./sum(nN_X,2),1:length(area2test),'v','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k')
for ar = 1 : length(area2test)
    text(0.05,ar,['Chi2=' num2str(chi2stat(ar)) ', p=' num2str(adj_p(ar))])
end    
set(gca,'Ytick',1:length(area2test),'YtickLabel',area2test,'FontSize',16)
        xlabel('Proportion of neurons')

%- 3 proportion across areas
nN = [nb_sig_reac(:,7) sum(nb_sig_reac,2) ];
x=0; 
clear pval_area chi2stat_area pairs
for ar1 = 1 : length(area2test)
    for ar2 = 1 : length(area2test)
        if ar1~=ar2 & ar1>ar2
           x=x+1
            [tbl,chi2stat_area(x),pval_area(x)] = chi2_fms(nN(ar1,1),nN(ar1,2),nN(ar2,1),nN(ar2,2));
            pairs(x,:)=[min([ar1 ar2]) max([ar1 ar2])]
        end
    end
end
[h, crit_p, adj_p]=fdr_bh(pval_area',.05);
a12l_p = adj_p(sum(pairs==4,2)==1)
a12l_chi = chi2stat_area(sum(pairs==4,2)==1)

%% Fig 5C - mixed selectivity at the Rew onset

params2take = [1 2];
mat2take = 'rew_diff';

%periods_reac = periods_all(2:3,:);
periods_mixed = {'REW' , 'Rew_onset' ,  [100 700]  };

time_chunk = bins_considered == find(ismember(times_evts,periods_mixed{1,2}));
time_considered = find(time_chunk & time>=periods_mixed{1,3}(1) & time<=periods_mixed{1,3}(2));

nb_sig_reac = [];nb_sig_reac_M = [];nb_sig_reac_X = [];
for ar = 1 : length(area2test)
    sig_reac = []; sig_reac_M = []; sig_reac_X = [];
    eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
    for p = 1 : length(params2take)
        dumm = squeeze(all_diff_sig(:,params2take(p),:))';
        sig_reac(:,p) = (sum(dumm(time_considered,takeit))~=0)';
        sig_reac_M(:,p) = (sum(dumm(time_considered,takeit & mk1))~=0)';
        sig_reac_X(:,p) = (sum(dumm(time_considered,takeit & ~mk1))~=0)';
    end
    nb_sig_reac(ar,:) = [sum(sig_reac(:,1)==1 & sig_reac(:,2)==0 ) ...
        sum(sig_reac(:,1)==0 & sig_reac(:,2)==1)  ...
        sum(sig_reac(:,1)==1 & sig_reac(:,2)==1 ) ...
        sum(sig_reac(:,1)==0 & sig_reac(:,2)==0 )];
    nb_sig_reac_M(ar,:) = [sum(sig_reac_M(:,1)==1 & sig_reac_M(:,2)==0 ) ...
        sum(sig_reac_M(:,1)==0 & sig_reac_M(:,2)==1)  ...
        sum(sig_reac_M(:,1)==1 & sig_reac_M(:,2)==1 ) ...
        sum(sig_reac_M(:,1)==0 & sig_reac_M(:,2)==0 )];
    nb_sig_reac_X(ar,:) = [sum(sig_reac_X(:,1)==1 & sig_reac_X(:,2)==0 ) ...
        sum(sig_reac_X(:,1)==0 & sig_reac_X(:,2)==1)  ...
        sum(sig_reac_X(:,1)==1 & sig_reac_X(:,2)==1 ) ...
        sum(sig_reac_X(:,1)==0 & sig_reac_X(:,2)==0 )];
end

     %   prop_sig=  nb_sig_reac./repmat(sum(nb_sig_reac')',1,size(nb_sig_reac,2)); %-  include the NS
        prop_sig= nb_sig_reac(:,1:end-1)./repmat(sum(nb_sig_reac(:,1:end-1)')',1,size(nb_sig_reac(:,1:end-1),2));
        prop_sig_M= nb_sig_reac_M(:,1:end-1)./repmat(sum(nb_sig_reac_M(:,1:end-1)')',1,size(nb_sig_reac_M(:,1:end-1),2));
        prop_sig_X= nb_sig_reac_X(:,1:end-1)./repmat(sum(nb_sig_reac_X(:,1:end-1)')',1,size(nb_sig_reac_X(:,1:end-1),2));

   figure;
   colorsBar = [cbrewer('qual', 'Set2', 8)];
   colorsBar = colorsBar([1 8 7],:);
    b = barh(prop_sig,'stacked'); hold on
    plot(prop_sig_M(:,1) ,0.2+(1:length(area2test)),'o','MarkerSize',10,'MarkerFaceColor',colorsBar(1,:),'MarkerEdgeColor','k')
    plot(prop_sig_M(:,1)+prop_sig_M(:,2) ,0.2+(1:length(area2test)),'o','MarkerSize',10,'MarkerFaceColor',colorsBar(3,:),'MarkerEdgeColor','k')
    plot(prop_sig_X(:,1)  ,-0.2+(1:length(area2test)),'v','MarkerSize',10,'MarkerFaceColor',colorsBar(1,:),'MarkerEdgeColor','k')
    plot(prop_sig_X(:,1)+prop_sig_X(:,2) ,-0.2+(1:length(area2test)),'v','MarkerSize',10,'MarkerFaceColor',colorsBar(3,:),'MarkerEdgeColor','k')

        box on;axis ij
        set(gca,'Ytick',1:length(area2test),'YtickLabel',area2test,'FontSize',16)
        xlabel('Proportion of neurons')
        xpos = [.05 .25 .55 .75];
        xlab = {'Rew' 'Rew+ChosenFlavor' 'ChosenFlavor' 'NS'};
        for i = 1 : size(prop_sig,2)
            b(i).FaceColor = colorsBar(i,:);
            text(xpos(i),.25,xlab{i},'Color',colorsBar(i,:),'FontSize',16)
        end

sig_reac = [];
        for p = 1 : length(params2take)
            dumm = squeeze(all_diff_sig(:,params2take(p),:))';
            sig_reac(:,p) = (sum(dumm(time_considered,:))~=0)';
        end
     sig = NaN(sum(keep),1);
     sig(   (sig_reac(:,1)==1 & sig_reac(:,2)==0) | (sig_reac(:,1)==0 & sig_reac(:,2)==1)   ,1)=false;
     sig(sig_reac(:,1)==1 & sig_reac(:,2)==1,1)=true;
 
     modeldata_su = table(sig(keep),all_units_area(keep),all_mk(keep),categorical(all_sess(keep)), 'VariableNames',{'sig' 'area' 'mk' 'sess'});
     modeldata_su(isnan(modeldata_su.sig),:)=[];
     modeldata_su.sig = (modeldata_su.sig)==1;

          models_form = {'sig ~ 1 + area  + (1|mk) + (1|sess)'};
            lme = fitglme(modeldata_su,models_form{1},'Distribution','Binomial'); 
            [~,wald_percNeurons,~,pval_adj_percNeurons]  = area_posthoc(lme,area2test,'y');
            disp(anova(lme));disp(pval_adj_percNeurons);disp(wald_percNeurons);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


 %- reward selective across areas
nN = [nb_sig_reac(:,1) sum(nb_sig_reac(:,1:end-1),2) ];
x=0; 
clear pval_area chi2stat_area pairs
for ar1 = 1 : length(area2test)
    for ar2 = 1 : length(area2test)
        if ar1~=ar2 & ar1>ar2
           x=x+1
            [tbl,chi2stat_area(x),pval_area(x)] = chi2_fms(nN(ar1,1),nN(ar1,2),nN(ar2,1),nN(ar2,2));
            pairs(x,:)=[min([ar1 ar2]) max([ar1 ar2])]
        end
    end
end
[h, crit_p, adj_p]=fdr_bh(pval_area',.05);
AI_p = adj_p(sum(pairs==8,2)==1)
AI_chi = chi2stat_area(sum(pairs==8,2)==1)
a12o_p = adj_p(sum(pairs==3,2)==1)
a12o_chi = chi2stat_area(sum(pairs==3,2)==1)


%% Main Figure ANOVA/LDA results (Figs 3-4)
figure;m = 0;

for p = 1 : length(params)
    m = m+1;
    
    % barplot percent sig
    for bins = 1:2 % : size(periods,1)
        time_chunk = bins_considered == find(ismember(times_evts,periods{bins,2}));
        time_considered = find(time_chunk & time>=periods{bins,3}(1) & time<=periods{bins,3}(2));
        for ar = 1 : length(area2test)
            eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
            dumm = squeeze(all_diff_sig(:,params(p),:))';
            
            bar_both(ar,:) = [sum(sum(dumm(time_considered,takeit))~=0) sum(takeit)];
            bar_M(ar,:) = [sum(sum(dumm(time_considered,takeit & mk1 ))~=0) sum(takeit & mk1)];
            bar_X(ar,:) = [sum(sum(dumm(time_considered,takeit & ~mk1 ))~=0) sum(takeit & ~mk1)];

            if bins == 2
                name_sig{p} = all_neuron_id(sum(dumm(time_considered,:))~=0 & keep',:);
            end
        end
        if bins == 1
            ref_sig = bar_both(:,1)./bar_both(:,2);
        end
        
        %- extract percent sig and stat
        if bins == 2
            sig = (sum(squeeze(all_diff_sig(keep,params(p),time_considered))')~=0)';
            
            modeldata_su = table(sig,all_units_area(keep),all_mk(keep),categorical(all_sess(keep)), 'VariableNames',{'sig' 'area' 'mk' 'sess'});
            models_form = {'sig ~ 1 + area  + (1|mk) + (1|sess)' ; 'sig ~ 1 + area  + (1|mk)'};
            % [lme,model_final] = model_comparison(modeldata_su,models_form,true);
            lme = fitglme(modeldata_su,models_form{1},'Distribution','Binomial'); model_final = models_form{1};

            [~,wald_percNeurons,~,pval_adj_percNeurons]  = area_posthoc(lme,area2test,'n');
            disp(['%%%%%%%%%% Percent of ' res_anova(1).(mat).predic{params(p)} ' neurons with model = ' model_final ' %%%%%%%%%%'])
            disp(anova(lme));disp(pval_adj_percNeurons);disp(wald_percNeurons);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        end
        
    end
    
    subplot(2,5,m)
    for ar = 1 : length(area2test)
        bar(ar, [bar_both(ar,1)./bar_both(ar,2)]*100,'FaceColor',colorsArea(ar,:));hold on
    end
    plot(ref_sig*100,'Color',[.6 .6 .6])
    plot([bar_M(:,1)./bar_M(:,2)]*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot([bar_X(:,1)./bar_X(:,2)]*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    ylim([0 80])
    xlim([0 length(area2test)+1])
    ylabel('Percent significant neurons')
    
    %- time resolved omega2
    m = m + 1;
    subplot(2,5,[m m+1])
    modeldata_o2=table(); %- for stats
    for ar = 1 : length(area2test)
        eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
        
        %- show perc sig
        % dumm = squeeze(all_diff_sig(:,params(p),:))';
        % perc = mean(dumm(:,takeit),2)*100;
        %- show Exp Variance in sig neurons
        time_chunk = bins_considered == find(ismember(times_evts,periods{2,2}));
        time_considered = find(time_chunk & time>=periods{2,3}(1) & time<=periods{2,3}(2));
        sig_anytime = ismember(bins_considered,timesel);
        
        sig = squeeze(all_diff_sig(:,params(p),:))';
        dumm = squeeze(all_diff_omega(:,params(p),:))';
        perc = nanmean(dumm(:,sum(sig(time_considered,:))>0 & takeit'),2);
        perc_sem = nanstd(dumm(:,sum(sig(time_considered,:))>0 & takeit')')' / sqrt(sum(takeit));

        om2(p,ar,1) = mean(mean(dumm(time_considered,sum(sig(time_considered,:))>0 & takeit'))   );
        om2(p,ar,2) = std(mean(dumm(time_considered,sum(sig(time_considered,:))>0 & takeit'))   )/sqrt(sum(takeit) );
        om2_m(p,ar) = mean(mean(dumm(time_considered,sum(sig(time_considered,:))>0 & takeit' & mk1'  ))   );
        om2_x(p,ar) = mean(mean(dumm(time_considered,sum(sig(time_considered,:))>0 & takeit' & ~mk1'  ))   );
       
        cd_curr = sum(sig(time_considered,:))>0 & takeit';
        modeldata_o2 = [modeldata_o2 ; table(double(mean(dumm(time_considered,cd_curr)))',all_units_area(cd_curr),all_mk(cd_curr),categorical(all_sess(cd_curr)), 'VariableNames',{'omega2' 'area' 'mk' 'sess'})];

        %- extract the Beta values
        sig1 = squeeze(all_diff_sig(:,params(1),:))';
        sig2 = squeeze(all_diff_sig(:,params(2),:))';
        dumm = squeeze(all_diff_betas(:,params(p),:))';
        dumm(~sig1 & ~sig2)=NaN;
        betas{p,ar} = nanmean(dumm(time_considered,takeit));
        betas_sigN{ar} = [sum(sig1(time_considered,takeit))>0 ; sum(sig2(time_considered,takeit))>0 ;...
                        sum(sig1(time_considered,takeit))>0 & sum(sig2(time_considered,takeit))>0 ];

        V = dumm(time_considered,sum(sig1(time_considered,:))>0 & sum(sig2(time_considered,:))>0 & takeit');
        betas_sig{p,ar} = mean(V);
        
        gaps = [1 find(diff(time)~=subsp)+1 length(time)];
        for t = 1 : length(timesel)
            time_considered = find(bins_considered == timesel(t) & time>=1000*timesel_sub(t,1) & time<=1000*timesel_sub(t,2));
            plot_me{t} = time_considered;
        end
        timestart = 1;
        xlab = [];
        for t = 1 : length(timesel)
            timeend = timestart+length(plot_me{t})-1;
            plot(timestart:timeend,perc(plot_me{t}),'Color',colorsArea(ar,:),'LineWidth',2); hold on
            plot(timestart:timeend,perc(plot_me{t})+perc_sem(plot_me{t}),'Color',colorsArea(ar,:),'LineWidth',.5); hold on
            plot(timestart:timeend,perc(plot_me{t})-perc_sem(plot_me{t}),'Color',colorsArea(ar,:),'LineWidth',.5); hold on
            line([timeend timeend],[0 80],'Color',[.6 .6 .6])
            timestart = timeend;
            xlab = [xlab time(plot_me{t})];
        end
        if strcmp(periods{2,1},'STIM')
            per = find(time(plot_me{1})>=periods{2,3}(1) & time(plot_me{1})<=periods{2,3}(2));
        else
            per = length(plot_me{1})+length(plot_me{2})+length(plot_me{3})+find(time(plot_me{4})>=periods{2,3}(1) & time(plot_me{4})<=periods{2,3}(2));
        end
        line([per(1) per(end)],[0 0],'Color','k','LineWidth',4)
           
        set(gca,'Xtick',1:step:length(xlab),'XtickLabel',xlab(1:step:end)/1000,'XtickLabelRotation',30,'FontSize',16)
        text(10,.1-(ar/200),[area2test{ar} ' - ' num2str(sum(takeit)) ],'Color',colorsArea(ar,:),'FontSize',16)
        
    end
    title(res_anova(1).(mat).predic{params(p)})
    ylabel('Mean Omega2')
    xlabel('Time (sec)')
    xlim([0 timeend+1])
    ylim([0 .1])
 
    %- stats
    models_form = {'omega2 ~ 1 + area  + (1|mk) + (1|sess)' ; 'omega2 ~ 1 + area  + (1|mk)'};
    %[lme,model_final] = model_comparison(modeldata_o2,models_form,false);
    lme = fitglme(modeldata_o2,models_form{1}); model_final = models_form{1};

    [~,wald_omegaNeurons,~,pval_adj_omegaNeurons]  = area_posthoc(lme,area2test,'n');
    disp(['%%%%%%%%%% Variance of ' res_anova(1).(mat).predic{params(p)} ' neurons with model = ' model_final ' %%%%%%%%%%'])
    disp(anova(lme));disp(pval_adj_omegaNeurons);disp(wald_omegaNeurons);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    if length(area2test)<=9
        %- decoding
        m = m +2;
        subplot(2,5,m);hold on
        
        dumm_areas = [];
        maxPerf = [];
        sessions_all=char();
        for ar = 1 : length(area2test)
            perf_all=[];
            nUnit=[];
            sessions=char();
            
            eval(['curr = LDA.res_LDA.' measures{p} '.' area2test{ar} ';']);
            for s = 1 : length(curr)
                nUnit(s) = length(curr(s).takeit);
                perf_all(s,:)=nanmean(curr(s).perf,2);
                sessions(s,:) = curr(s).lda_sess(1:7);
            end
            
            %- decoding only done on bin 2 = so time to take [200-800]
            %t_lim = length(time_considered);
            
            nUnit_all{p,ar} = nUnit;
            %  maxPerf = [maxPerf ; [nanmean(perf_all(:,end-t_lim:end)')' , ar*ones(length(nUnit),1) , nUnit']];
            maxPerf = [maxPerf ; [perf_all , ar*ones(length(nUnit),1) , nUnit']];
            sessions_all = [sessions_all ; sessions];
            mk = sessions_all(:,1);
            mk = ismember(mk,'M')+1;
            dumm_areas = [dumm_areas ; repmat(area2test(ar),length(nUnit),1)];
        end
                
        %- plot decoding perf
        mm = 0;
        for i  = 1 : length( area2test)
            X = maxPerf(maxPerf(:,2)==i,1);
            mk_sub = mk(maxPerf(:,2)==i);
            yl=i+mm;
            wdth = .5;
         %   boxplot_ind_mk(X,yl,wdth,[.8 .8 .8 ; .65 .65 .65 ; colorsArea(i,:)],mk_sub)
               boxplot_ind(X,yl,wdth,[colorsArea_sub(i,:) ; colorsArea(i,:)]);
               nbsess(i,:) = [sum(mk_sub==2) sum(mk_sub==1)];

        end
        set(gca,'view',[90 -90],'color','none','FontSize',16);
        set(gca,'YTick',1:length(area2test),'YTickLabel',area2test,'YTickLabelRotation',30)
        ylim([0 length(area2test)+1]);
        title(measures{p})
        xlabel('Decoding probability')
        
        %- significance
        modeldata_pop = table(maxPerf(:,1),dumm_areas,maxPerf(:,3),sessions_all,sessions_all(:,1),'VariableNames',{'perf' 'area' 'nb' 'sess' 'mk'});
        models_form = {'perf ~ 1 + area  + (1|mk) + (1|sess)' ; 'perf ~ 1 + area  + (1|mk)'};
        % [lme,model_final] = model_comparison(modeldata_pop,models_form,false);
        lme = fitglme(modeldata_pop,models_form{1}); model_final = models_form{1};

        [pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'n');
        disp(['%%%%%%%%%% Decoding perf with model = ' model_final ' %%%%%%%%%%'])
        disp(anova(lme));disp(pval_adj);disp(wald);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        %- plot the significance
        pval2=[pval_adj , zeros(length(area2test),1)]+ [pval_adj' ; zeros(1,length(area2test))];
        for ar1  = 1 : length( area2test)
            Xmax = max(maxPerf(maxPerf(:,2)==ar1,1));
            Xmin = min(maxPerf(maxPerf(:,2)==ar1,1));
            Xmean = mean(maxPerf(maxPerf(:,2)==ar1,1));
            updown=[1.5 2];
            for ar2 = 1 : length(area2test)
                Xmean2 = mean(maxPerf(maxPerf(:,2)==ar2,1));
                if ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean<Xmean2
                    text(Xmax+(updown(1)/120),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
                    updown(1) = updown(1) + 1;
                elseif ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean>Xmean2
                    text(Xmin-(updown(2)/120),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
                    updown(2) = updown(2) + 1;
                end
            end
        end
        
        %- sig decoding proportion
        m = m + 1;
        subplot(2,5,m);hold on
        
        dumm_perf=[];dumm_mk =[];dumm_sess=[];dumm_area=[];
        for ar = 1 : length(area2test)
            eval(['curr = LDA.res_LDA.' measures{p} '.' area2test{ar} ';']);
            perf = [curr(:).perf_pval];
            bar_decod(ar,:) = [sum(perf<0.05) length(perf<0.05) ];
            
            sess = cat(1,curr(:).lda_sess);
            mk = sess(:,1);
            
            bar_decodM(ar,:) = [sum(perf(mk=='M')<0.05) length(perf(mk=='M')) ];
            bar_decodX(ar,:) = [sum(perf(mk=='X')<0.05) length(perf(mk=='X')) ];
            dumm_perf = [dumm_perf ; perf'<0.05 ];
            dumm_mk = [dumm_mk ; mk ];
            dumm_sess = [dumm_sess ; sess ];
            dumm_area = [dumm_area ; repmat(area2test(ar),size(mk)) ];
        end
        for ar = 1 : length(area2test)
            bar(ar, [bar_decod(ar,1)./bar_decod(ar,2)]*100,'FaceColor',colorsArea(ar,:));hold on
        end
        plot([bar_decodM(:,1)./bar_decodM(:,2)]*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
        plot([bar_decodX(:,1)./bar_decodX(:,2)]*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
        set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
        ylim([0 100])
        xlim([0 length(area2test)+1])
        ylabel('Percent significant decoding')
        
        %- significance percent sig decoding
            
        modeldata_pop = table(dumm_perf,dumm_area,dumm_mk,dumm_sess, 'VariableNames',{'sig' 'area' 'mk' 'sess'});
        models_form = {'sig ~ 1 + area  + (1|mk) + (1|sess)' ; 'sig ~ 1 + area  + (1|mk)'};
        % [lme,model_final] = model_comparison(modeldata_pop,models_form,true);
        lme = fitglme(modeldata_pop,models_form{1},'Distribution','Binomial'); model_final = models_form{1};

        [pval,wald_popsig_adj,thr_corr,pval_popsig_adj] = area_posthoc(lme,area2test,'n');
        disp(['%%%%%%%%%% Percent sig decoding with model = ' model_final ' %%%%%%%%%%'])
        disp(anova(lme));disp(pval_popsig_adj);disp(wald_popsig_adj);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        %- Build Table s2
        tab_s2_part1 =[];ar_s2 =[];
        for ar1 = 1 : length(area2test)-1
            for ar2 = 1 : length(area2test)
                if ar2>ar1
                    n=1;pp=0;
                    while pp == 0
                        n = n + 1;
                        pp = round(pval_adj_percNeurons(ar2,ar1),n);
                    end
                    pp = round(pval_adj_percNeurons(ar2,ar1),n+1);
                   
                    tab_s2_part1 = [tab_s2_part1;{num2str(round(wald_percNeurons(ar2,ar1),2)) num2str(pp)}];
                    ar_s2 = [ar_s2 ;{area2test{ar1} area2test{ar2}}];
                end
            end
        end
        tab_s2_part2 =[];ar_s2 =[];
        for ar1 = 1 : length(area2test)-1
            for ar2 = 1 : length(area2test)
                if ar2>ar1
                    n=1;pp=0;
                    while pp == 0
                        n = n + 1;
                        pp = round(pval_adj(ar2,ar1),n);
                    end
                    pp = round(pval_adj(ar2,ar1),n+1);
                   
                    tab_s2_part2 = [tab_s2_part2;{num2str(round(wald(ar2,ar1),2)) num2str(pp)}];
                    ar_s2 = [ar_s2 ;{area2test{ar1} area2test{ar2}}];
                end
            end
        end
        tab_s2_part3 =[];ar_s3 =[];
        for ar1 = 1 : length(area2test)-1
            for ar2 = 1 : length(area2test)
                if ar2>ar1
                    n=1;pp=0;
                    while pp == 0
                        n = n + 1;
                        pp = round(pval_adj_omegaNeurons(ar2,ar1),n);
                    end
                    pp = round(pval_adj_omegaNeurons(ar2,ar1),n+1);
                   
                    tab_s2_part3 = [tab_s2_part3;{num2str(round(wald_omegaNeurons(ar2,ar1),2)) num2str(pp)}];
                    ar_s3 = [ar_s3 ;{area2test{ar1} area2test{ar2}}];
                end
            end
        end

        if p == 1
            tab_s2 = [tab_s2_part1 tab_s2_part3 tab_s2_part2];
        else
            tab_s2 = [tab_s2 tab_s2_part1 tab_s2_part3 tab_s2_part2];
        end
        
    else
        m = m + 3;

    end
    
end
set(gcf, 'Color', [1 1 1]);

% save('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\POTT_sigUnits_name.mat','name_sig')

grpstats(modeldata_pop,{'area' 'mk'},"numel")

figure;
for m = 1 : 2
    subplot(1,2,m)
    for ar = 1 : length(area2test)
        bar(ar, squeeze(om2(m,ar,1)) ,'FaceColor',colorsArea(ar,:));hold on
        line([ar ar], [squeeze(om2(m,ar,1))-squeeze(om2(m,ar,2)) squeeze(om2(m,ar,1))+squeeze(om2(m,ar,2))] ,'Color',colorsArea(ar,:));hold on
        plot(ar,om2_m(m,ar),'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
        plot(ar,om2_x(m,ar),'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    end
    ylim([0 0.07])    
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    xlim([0 length(area2test)+1])
    ylabel('Mean Omega2')

end
