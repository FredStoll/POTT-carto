%- decrease nb trial used to 5 instead of 10! otherwise nothing for Morbier

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
list = dir([path2go 'POOL_2*.mat']);
params = load([path2go 'M021519a_SPKpool.mat']);

param.predic = {'I_chosenproba_by_juice' 'I_juice_by_chosenproba'}; %- select param you want to test..
param.minNeurons = [25 50 75 100 150 200 300 400 500 600 700 800 900 1000]; %- min number of neuron to run
param.nComp = 20; %- nb of component to keep for PCA
param.Repetition = 200; %- number of times to run the decoder (100 or 200 ideally)
param.minTrCond = 5; %- minimum number of trials per condition to run
param.overwrite = true;
param.thresh = true;
param.thr = 1; % in Hz, used if param.thresh is 'y'
param.rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
param.bins4decoding = [2] ; %- 6 for rew %- perform decoding on subset on bins (stim period here)
param.timebin = [100 700] ;% 0 600 rew %- perform decoding on subset on bins (stim period here)
param.normalize_fr = false;
param.tr_perm = true;

disp(['Computing LDA on ' num2str(length(list)) ' sessions'])
areas = utils_POTT_areas;

area2test = {'12r' '12m' '12o' '12l' 'a11ml' '13l' '13m' 'LAI'};

util.minmaxnorm = @(data) (data-min(data))/(max(data)-min(data));
util.minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

un = find(list(1).name=='_');

%- check if matrix already exists for each predictor
if exist([path2go 'res_LDAcrosspop_kfold_final_' list(1).name(un(1)+1:un(2)-1) '.mat'])==2
    prev = load([path2go 'res_LDAcrosspop_kfold_final_' list(1).name(un(1)+1:un(2)-1) '.mat']);
    done = [];
    for pr = 1 : length(param.predic)
        if isfield(prev.res_LDAcrosspop,param.predic{pr})
            done(pr) = true;
        end
    end
    if ~param.overwrite
        param.predic(find(done==true))=[];
    end
    res_LDApop = prev.res_LDAcrosspop;
end
if param.overwrite
    clear res_LDAcrosspop
end

if ~isempty(param.predic) %- if there is some predictors not done yet

    for pr = 1 : length(param.predic)
        %nSess = zeros(length(area2test),1);
        for s = 1 : length(list)
            clearvars -except path2go list param params areas area2test util s u pr res_LDAcrosspop nSess un
            name = list(s).name;

            %- load spike data and histology info
            load([path2go name]); disp(s);


            %- create time vector and remove overlapping windows/border events
            bins2remove = (param.rmv/params.subsp)-1;

            n_evt = length(params.times_evts);
            lin=size(neurons_info,1)*n_evt;
            nTimes = [0 ; (sum(abs(params.times(:,:)),2)/params.subsp)];
            col = sum(nTimes);
            bins = NaN(1,col);
            time = NaN(1,col);

            for b = 1 : length(params.times_evts)
                time(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) =  (params.times(b,1):params.subsp:params.times(b,2)-(params.subsp/2));
                bins(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) = b;
                bins([sum(nTimes(1:b))+1:sum(nTimes(1:b))+1+bins2remove   sum(nTimes(1:b+1))-bins2remove:sum(nTimes(1:b+1))    ]) = NaN; %- remove 100 ms each side (avoid overlaps between bins and smoothing problems)
            end

            %- extract the data for the task 2 consider
            if strcmp(param.predic{pr}(1),'I') %- instrumental
                SPK_data = SPK_INS;
                BEHAV = BEHAV_INS;
            elseif strcmp(param.predic{pr}(1),'P') %- pav
                SPK_data = SPK_PAV;
                BEHAV = BEHAV_PAV;
            end

            %- normalize FR + cut to only bin 4 decoding
            if param.normalize_fr %- normalize FR using min-max in considered bins only (all trials except border of events)
                SPK_data_norm = SPK_data;
                units = unique(Tr_Clust_data(:,2));
                for u = 1 : length(units)
                    temp = SPK_data(Tr_Clust_data(:,2)==units(u),~isnan(bins)); %- could normalize only on the bin2use, but more likely to get some 0 FR?? so NaNs!
                    SPK_data_norm(Tr_Clust_data(:,2)==units(u),~isnan(bins)) = reshape(util.minnorm(temp(:)),size(temp,1),size(temp,2));
                end

                SPK_data_cut = SPK_data_norm;
                SPK_data_cut(:,~ismember(bins,param.bins4decoding))=[];
            else
                SPK_data_cut = SPK_data;
                SPK_data_cut(:,~ismember(bins,param.bins4decoding))=[];
            end

            time(~ismember(bins,param.bins4decoding)) = [];

            %- reject low FR neurons, if needed
            if param.thresh  %- reject neurons with FR too low
                reject = [];
                for n = 1 : length(unique(BEHAV(1,:)))
                    if mean(mean(SPK_data_cut(BEHAV(1,:)==n,:)))<param.thr  %- put NaNs when FR is below thr
                        reject = [reject , n];
                    end
                end

                %- remove them from main matrices
                SPK_data_cut( ismember(BEHAV(1,:),reject) ,:) = [];
                BEHAV(:, ismember(BEHAV(1,:),reject) ) = [];
            end

            %- create data and factor matrix for the decoding
            data = SPK_data_cut;
            unit_nb = unique(BEHAV(1,:));

            TrialType_header = ['Unit_Nb' 'Session_Nb' params.TrialType_header];

            if sum(ismember(TrialType_header,param.predic{pr})) ~= 0
                dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' param.predic{pr} }),:)';
            else %- for the proba/juice interaction
                dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' 'I_chosenproba' 'I_chosenjuice' }),:)';
                dumm = [dumm(:,1) dumm(:,2)+(100*dumm(:,3))];
            end

            factor = dumm(:,[2 1]);

            %- trials to consider for instrum
            if strcmp(param.predic{pr}(1),'I')
                when = @(arg,cd) BEHAV(strcmp(TrialType_header,arg),:)==cd ;
                diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
                %    diff_juice_all = repmat(diff_juice',length(unique(unit_nb)),1);
            end

            %- only keep a subset of proba
            remove = [];
            if strcmp(param.predic{pr},'P_proba')
                keep = [10 30 50 70 90 ];
                remove = ~ismember(dumm(:,2),keep);
            elseif strcmp(param.predic{pr},'I_chosenproba')
                keep = [30 50 70 90 ];
                remove = ~ismember(dumm(:,2),keep) | ~diff_juice';  % was not removing same juice trials before
            elseif strcmp(param.predic{pr},'I_unchosenproba')
                keep = [10 30 50 70];
                remove = ~ismember(dumm(:,2),keep) | ~diff_juice';  % was not removing same juice trials before
            elseif strcmp(param.predic{pr},'I_chosenjuice')
                keep = [1 2];
                remove = ~diff_juice';  % was not removing same juice trials before
            elseif strcmp(param.predic{pr},'I_chosenside')
                keep = [1 2];
                remove = ~diff_juice';  % was not removing same juice trials before
            elseif strcmp(param.predic{pr},'I_chosenproba_juice') | strcmp(param.predic{pr},'I_chosenproba_by_juice') | strcmp(param.predic{pr},'I_juice_by_chosenproba')
                keep = [130 150 170 190 230 250 270 290];
                remove = ~ismember(dumm(:,2),keep) | ~diff_juice';
            else
                keep = unique(dumm(:,2));
            end
            if ~isempty(remove)
                data(remove,:) = [];
                dumm(remove,:)=[];
                factor(remove,:)=[];
            end

            %- check if enough trial and remove neurons without enough
            for p = 1 : length(keep)
                for n = 1 : length(unit_nb)
                    nTr(n,p) = sum(dumm(:,2)==keep(p) & dumm(:,1)==unit_nb(n));
                end
            end

            remove = sum(nTr>param.minTrCond,2)~= length(keep);
            unit_nb(remove)=[];

            %- extract only 1 bin during stim onset
            data = mean(data(:,time>=param.timebin(1) & time<=param.timebin(2)),2);

            for u = 1 : length(param.minNeurons)
                if length(unit_nb)>=param.minNeurons(u)
                    for p = 1 : param.Repetition
                        disp(['Area ' num2str(s) '/' num2str(length(list)) ' - ' num2str(param.minNeurons(u)) ' units - ' param.predic{pr} ' - perm ' num2str(p) '/' num2str(param.Repetition)])

                        if ~isempty(param.minNeurons(u))
                            if param.minNeurons(u)<=length(unit_nb)
                                rd = randperm(length(unit_nb));
                                unit_nb_sub = unit_nb(rd);
                                unit_nb_sub = unit_nb_sub(1:param.minNeurons(u));
                            end
                        end

                        data_sub =  data( ismember(factor(:,2),unit_nb_sub),:);
                        factor_sub =  factor( ismember(factor(:,2),unit_nb_sub),:);

                        if strcmp(param.predic{pr},'I_chosenproba_by_juice') %- in that case, run 2 decoders, 1 per juice!

                            %-  juice 1 trials
                            data_sub_j1 = data_sub;
                            factor_sub_j1 = factor_sub;
                            j1 = [130 150 170 190];
                            remove = ~ismember(factor_sub_j1(:,1),j1);
                            data_sub_j1(remove,:) = [];
                            factor_sub_j1(remove,:)=[];

                            [XX,YY,param_LDA] = pop_init_noRdmTrials_FORCEDMIN(data_sub_j1,factor_sub_j1,1,'perm',1,'minTr',param.minTrCond,'pop',false); % ,'window',[-.5 1.5]

                            %-  juice 2 trials
                            data_sub_j2 = data_sub;
                            factor_sub_j2 = factor_sub;
                            j2 = [230 250 270 290];
                            remove = ~ismember(factor_sub_j2(:,1),j2);
                            data_sub_j2(remove,:) = [];
                            factor_sub_j2(remove,:)=[];

                            [XX2,YY2,param_LDA_2] = pop_init_noRdmTrials_FORCEDMIN(data_sub_j2,factor_sub_j2,1,'perm',1,'minTr',param.minTrCond,'pop',false); % ,'window',[-.5 1.5]

                            [perf(:,p),~] = pca_lda_kfold_crosspop(XX(1),YY(1),XX2(1),YY2(1),param);

                            %- randomize the 2 groups so that we can do
                            %- decoding with same nb of trials but no
                            %- cross-condition anymore
                            XX_mix{1} = [];
                            XX2_mix{1} = [];
                            bothXX = [XX{1} ; XX2{1}];
                            both = [YY{1} ; YY2{1}];
                            both_grp = mod(both,100);
                            cds = unique(both_grp);
                            for cd = 1 : length(cds)
                                tr_dumm = find(both_grp==cds(cd));
                                training_testing = randperm(size(tr_dumm,1));
                                XX_mix{1} = [XX_mix{1} ; bothXX(tr_dumm(training_testing(1:size(tr_dumm,1)/2)),:)];
                                XX2_mix{1} = [XX2_mix{1} ; bothXX(tr_dumm(training_testing((size(tr_dumm,1)/2)+1:end)),:)];
                            end
                            [perf_raw(:,p),~] = pca_lda_kfold_crosspop(XX_mix,YY,XX2_mix,YY2,param);

                            yy_perm = randperm(length(YY{1}));
                            for tt = 1 : length(YY)
                                YYperm{tt} = YY{1}(yy_perm);
                            end
                            [perf_perm(:,p),~] = pca_lda_kfold_crosspop(XX,YYperm,XX2(1),YY2(1),param);

                        elseif strcmp(param.predic{pr},'I_juice_by_chosenproba') %- in that case, run 2 decoders, 1 per juice!

                            [XX,YY,param_LDA] = pop_init_noRdmTrials_FORCEDMIN(data_sub,factor_sub,1,'perm',1,'minTr',param.minTrCond,'pop',false); % ,'window',[-.5 1.5]

                            pairs ={[130 230 ; 150 250] ; [130 230 ; 170 270] ; [130 230 ; 190 290] ; ...
                                [150 250 ; 170 270] ; [150 250 ; 190 290] ; ...
                                [170 270 ; 190 290]};

                            x = -1;
                            for co = 1 : length(pairs)
                                x = x + 2;
                                YY1 = ismember(YY{1},pairs{co}(1,:));
                                YY2 = ismember(YY{1},pairs{co}(2,:));
                                [temp_perf(x:x+1),~] = pca_lda_kfold_crosspop_juice({XX{1}(YY1,:)},{YY{1}(YY1)},{XX{1}(YY2,:)},{YY{1}(YY2)},param);

                                %- randomize the 2 groups so that we can do
                                %- decoding with same nb of trials but no
                                %- cross-condition anymore
                                XX_mix{1} = [];
                                XX2_mix{1} = [];
                                bothXX = [XX{1}(YY1,:) ; XX{1}(YY2,:)];
                                both = [YY{1}(YY1) ; YY{1}(YY2)];
                                both_grp = floor(both/100);
                                cds = unique(both_grp);
                                for cd = 1 : length(cds)
                                    tr_dumm = find(both_grp==cds(cd));
                                    training_testing = randperm(size(tr_dumm,1));
                                    XX_mix{1} = [XX_mix{1} ; bothXX(tr_dumm(training_testing(1:size(tr_dumm,1)/2)),:)];
                                    XX2_mix{1} = [XX2_mix{1} ; bothXX(tr_dumm(training_testing((size(tr_dumm,1)/2)+1:end)),:)];
                                end
                                [temp_perf_raw(x:x+1),~] = pca_lda_kfold_crosspop_juice(XX_mix,{YY{1}(YY1)},XX2_mix,{YY{1}(YY2)},param);


                                yy_perm = randperm(length(YY{1}));
                                for tt = 1 : length(YY)
                                    YYperm{tt} = YY{1}(yy_perm);
                                end
                                YY1 = ismember(YYperm{1},pairs{co}(1,:));
                                YY2 = ismember(YYperm{1},pairs{co}(2,:));
                                [temp_perf_perm(x:x+1),~] = pca_lda_kfold_crosspop_juice({XX{1}(YY1,:)},{YYperm{1}(YY1)},{XX{1}(YY2,:)},{YYperm{1}(YY2)},param);


                            end

                            perf(1,p) = nanmean(temp_perf);
                            perf_raw(1,p) = nanmean(temp_perf_raw);
                            perf_perm(1,p) = nanmean(temp_perf_perm);
                        end
                    end

                    %- extract param of interest to save
                    %  nSess(ar,1) = nSess(ar,1) + 1 ;
                    temp.lda_sess = name;
                    temp.perf = perf;
                    temp.perf_raw = perf_raw;
                    temp.nbUnit = param.minNeurons(u);
                    if param.tr_perm
                        temp.perf_perm = perf_perm;
                    end
                    unde = find(name=='_',1,'last');
                    area_name = ['area_' name(unde+1:end-4)];
                    eval(['res_LDAcrosspop.' param.predic{pr} '.' area_name '(u) = temp;'])
                end
            end
        end
        eval(['res_LDAcrosspop.' param.predic{pr} '.param = param;'])
    end

    save([path2go 'res_LDAcrosspop_kfold_final_' list(1).name(un(1)+1:un(2)-1) '.mat'],'res_LDAcrosspop','list','area2test','path2go')
end

%% POST HOC - FIGURE

clear

nbNeur = 100;

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
load([path2go 'res_LDAcrosspop_kfold_final_289sessions.mat'])
measures = {'I_chosenproba_by_juice' 'I_juice_by_chosenproba'}

chance_l = [.25 .5];
evtCol = {'Greens' 'Greens'};

colorsArea = cbrewer('qual', 'Paired', 12);
if length(area2test)==5
    colorsArea = colorsArea([2 4 5 6 12],:);
else
    colorsArea = colorsArea([1:7 12 8:11],:);
end

meanPerf={};
zval_diff0=[]; pval_diff0=[];
for m = 1 : length(measures)
    x=0;
    for ar = 1: length(area2test)
        eval(['temp = res_LDAcrosspop.' measures{m} '.area_' area2test{ar} ';'])
        for n = 1 : length(temp)
            if temp(n).nbUnit<=500
                %  temp(n).perf = temp(n).perf_raw;
                x = x+1;
                %pval = sum(temp(n).perf_perm > nanmean(temp(n).perf)) / length(temp(n).perf_perm);
                sem = nanstd(temp(n).perf);%/sqrt(length(temp(n).perf)));
                sem_raw = nanstd(temp(n).perf_raw);%/sqrt(length(temp(n).perf_raw)));
                sem_perm = nanstd(temp(n).perf_perm);%/sqrt(length(temp(n).perf_perm)));
                meanPerf{m}(x,:) = [nanmean(temp(n).perf) sem nanmean(temp(n).perf_raw) sem_raw nanmean(temp(n).perf_perm)  sem_perm prctile((temp(n).perf_raw-temp(n).perf)./temp(n).perf_raw,[50 25 75]) temp(n).nbUnit ar ];
            end
            if temp(n).nbUnit==nbNeur
                y_diff = (temp(n).perf_raw-temp(n).perf)./temp(n).perf_raw;
                y_diff(isnan(y_diff))=[];
                [pval_diff0(m,ar),h,stats] = signrank(y_diff);
                zval_diff0(m,ar) = stats.zval;
            end
        end
    end
    meanPerf_sub{m} = meanPerf{m}(meanPerf{m}(:,end-1)==nbNeur,:);

end
[h, crit_p, adj_p]=fdr_bh(pval_diff0,.05);

xax = repmat(meanPerf_sub{m}(:,end)',2,1)+[-.15  .15]';
figure;
for m = 1 : length(measures)
    subplot(1,2,m)
    pp = plot(xax,[meanPerf_sub{m}(:,3) meanPerf_sub{m}(:,1) ]','.-','MarkerSize',35,'LineWidth',2);hold on
    for ar = 1 : length(pp)
        set(pp(ar),'Color',colorsArea(meanPerf_sub{m}(ar,end),:))
        line([xax(1,ar) xax(1,ar)],[meanPerf_sub{m}(ar,3)-meanPerf_sub{m}(ar,4) meanPerf_sub{m}(ar,3)+meanPerf_sub{m}(ar,4)],'Color',colorsArea(meanPerf_sub{m}(ar,end),:))
        line([xax(2,ar) xax(2,ar)],[meanPerf_sub{m}(ar,1)-meanPerf_sub{m}(ar,2) meanPerf_sub{m}(ar,1)+meanPerf_sub{m}(ar,2)],'Color',colorsArea(meanPerf_sub{m}(ar,end),:))
    end
    ylim([chance_l(m)-0.05 1])
    xlim([0 length(area2test)+1])
    set(gca,'Xtick',1:length(area2test),'XTickLabel',area2test,'Fontsize',18)

end

xax = repmat(meanPerf_sub{m}(:,end)',2,1)+[-.15  .15]';
figure;
for m = 1 : length(measures)
    subplot(2,1,m)
    for ar = 1 : length(area2test)
        plot(100*meanPerf_sub{m}(ar,7),ar,'.','MarkerSize',35,'LineWidth',2,'Color',colorsArea(meanPerf_sub{m}(ar,end),:));hold on
        line(100*[meanPerf_sub{m}(ar,8) meanPerf_sub{m}(ar,9)],[ar ar],'Color',colorsArea(meanPerf_sub{m}(ar,end),:))
    end
    ylim([0 length(area2test)+1])
    set(gca,'Ytick',1:length(area2test),'YTickLabel',area2test,'Xtick',[-40:20:40],'XTickLabel',[-40:20:40],'Fontsize',18)
    xlim(100*[-.4 .4])
    axis ij
end

%- Figure 7B-C
all_decoding_perf = 0;
for  m = 1  : length(measures)
    % fig = figure('Position',[1051 51 375 918]);
    fig = figure('Position',[1125 51 301 918]);
    subplot(3,1,[1 2]);
    for ar = 1: length(area2test)
        dumm = meanPerf{m}(meanPerf{m}(:,end)==ar,[1 end-1 2]);
        [param_a(m,ar),param_b(m,ar),Yh] = fitDecodPerf(dumm,chance_l(m));
        plot(dumm(:,2),dumm(:,1),'.','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15,'MarkerFaceColor',colorsArea(ar,:),'MarkerEdgeColor',colorsArea(ar,:));hold on
        plot(dumm(:,2),Yh+chance_l(m),'-','Color',colorsArea(ar,:),'LineWidth',1,'MarkerSize',15);hold on

        dumm = meanPerf{m}(meanPerf{m}(:,end)==ar,[3 end-1 4]);
        [param_a(m,ar),param_b(m,ar),Yh] = fitDecodPerf(dumm,chance_l(m));
        plot(dumm(:,2),Yh+chance_l(m),'--','Color',colorsArea(ar,:),'LineWidth',1,'MarkerSize',15);hold on

        keep4legend(ar,1) = dumm(1,1);

        % text(825,1-(ar/60),area2test{ar},'Color',colorsArea(ar,:),'FontWeight','bold','FontSize',16)
        if chance_l(m)<0.5
            ylim([0.2 1])
        else
            ylim([0.4 0.9])
        end
        xlim ([-150 600])
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
end
