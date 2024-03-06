%%  PCA + subspace PCA - Figures 7/9
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.09

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
list = dir([path2go 'POOL_2*.mat']);
params = load([path2go 'M021519a_SPKpool.mat']);

param.predic = {'I_chosenproba_by_juice' } % 'I_chosenside'}; %- select param you want to test..
%param.predic = {'P_juice' 'P_proba' 'P_side'}; %- select param you want to test..
param.minNeurons = [50 100 150 200]; %- min number of neuron to run
param.nComp = 20; %- nb of component to keep for PCA
param.Repetition = 200; %- number of times to run the decoder (100 or 200 ideally)
param.minTrCond = 5; %- minimum number of trials per condition to run
param.overwrite = false;
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
util.meannorm  = @(data) (data-mean(data))/(max(data)-min(data));
util.minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

for pr = 1 : length(param.predic)
    %nSess = zeros(length(area2test),1);
    for s = 1 : length(list)
        clearvars -except path2go list param params areas area2test util s u pr res_PCA nSess un ALLdata
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

        %  predic_all = repmat(dumm(:,2),length(unique(unit_nb)),1);
        %  factor = [dumm(:,2) , unit_nb];

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

        conditions = unique(factor(:,1));
        for n = 1 : length(unit_nb)
            data_sub = data( ismember(factor(:,2),unit_nb(n)),:);
            factor_sub = factor(ismember(factor(:,2),unit_nb(n)),1);

            data_norm = (data_sub-mean(data_sub(:)))/std(data_sub(:));
            data_centered = data_norm-mean(data_norm);

            for f = 1 : length(conditions)
                ALLdata{s}(n,f,:) = mean(data_centered(factor_sub==conditions(f),:));
            end
        end
    end
end

%- plot PCA for Fig 7
figure;
for ar = 1 : length(area2test)
    neuron_selection =  randperm(size(ALLdata{ar},1),100);

    ALL = [squeeze(ALLdata{ar}(neuron_selection,1,:)) , ...
        squeeze(ALLdata{ar}(neuron_selection,2,:)) , ...
        squeeze(ALLdata{ar}(neuron_selection,3,:)) , ...
        squeeze(ALLdata{ar}(neuron_selection,4,:)) , ...
        squeeze(ALLdata{ar}(neuron_selection,5,:)) , ...
        squeeze(ALLdata{ar}(neuron_selection,6,:)) , ...
        squeeze(ALLdata{ar}(neuron_selection,7,:)) , ...
        squeeze(ALLdata{ar}(neuron_selection,8,:)) ]';

    rejected = find(sum(isnan(ALL))>0);
    %disp(['Removing ' num2str(length(rejected)) ' neurons (FR < ' num2str(thr) ' Hz)'])
    ALL(:,sum(isnan(ALL))>0) = [];
    ALL(:,sum(~isfinite(ALL))>0) = [];

    [coeff,score,~,~,explained,~] = pca(ALL,'NumComponents',3);

    score1 = score(1:length(time),:);
    score2 = score(length(time)+1:2*length(time),:);
    score3 = score(2*length(time)+1:3*length(time),:);
    score4 = score(3*length(time)+1:4*length(time),:);
    score5 = score(4*length(time)+1:5*length(time),:);
    score6= score(5*length(time)+1:6*length(time),:);
    score7 = score(6*length(time)+1:7*length(time),:);
    score8 = score(7*length(time)+1:end,:);

    %- plot the PCA result
    %colors = cbrewer('div', 'RdBu', 20);
    % colors = flipud(colors);
    colors = cbrewer('div', 'PiYG', 8);
    colors = colors([1:4 8 7 6 5],:);

    subplot(2,4,ar)
    hold on
    minT=0;
    maxT=700;
    for ii = 1 : 8
        eval(['score_curr = score' num2str(ii) ';'])
        h = plot3(score_curr(time>=minT & time<=maxT,1),score_curr(time>=minT & time<=maxT,2),score_curr(time>=minT & time<=maxT,3))%,'.-','MarkerSize',15);
        xlabel(['pc1 (' num2str(explained(1)) '%)']) ;ylabel(['pc2 (' num2str(explained(2)) '%)']); zlabel(['pc3 (' num2str(explained(3)) '%)'])
        set(h, {'color'}, {[colors(ii,:)]});
        loc = double(score_curr(time>minT & time<maxT,1:3))-0.001;
        if ii == 1
            text(loc(1,1),loc(1,2),loc(1,3),'Stim')
        end
        h = plot3(loc(1,1)+0.001,loc(1,2)+0.001,loc(1,3)+0.001,'.','MarkerSize',15);
        set(h, {'color'}, {[colors(ii,:)]});
    end
    title(area2test{ar})

    perc_expl(ar)=sum(explained(1:3));
end


%% pre-processing for subspace PCA (Fig 9)
clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
list = dir([path2go 'POOL_2*.mat']);
params = load([path2go 'M021519a_SPKpool.mat']);

param.predic = {'I_proba_rew' } % 'I_chosenside'}; %- select param you want to test..
param.minTrCond = 5; %- minimum number of trials per condition to run
param.overwrite = false;
param.thresh = true;
param.thr = 1; % in Hz, used if param.thresh is 'y'
param.rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
param.bins4decoding = [2 4 5 6] ; %- 6 for rew %- perform decoding on subset on bins (stim period here)
param.normalize_fr = false;
param.tr_perm = true;

disp(['Computing LDA on ' num2str(length(list)) ' sessions'])
areas = utils_POTT_areas;

area2test = {'12r' '12m' '12o' '12l' 'a11ml' '13l' '13m' 'LAI'};

util.minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

un = find(list(1).name=='_');
%- check if matrix already exists (big ~3Gb file, don't put on the cloud!)
if exist('C:\Users\Fred\Documents\POTTc_PCA_subspaces_perm.mat')==2 & ~param.overwrite
    param.predic = [];
end

if ~isempty(param.predic) %- if there is some predictors not done yet

    for pr = 1 : length(param.predic)
        %nSess = zeros(length(area2test),1);
        for s = 1 : length(list)
            clearvars -except path2go list param params areas area2test util s u pr res_PCA nSess un ALLdata REWdata PERMdata PBdata
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
                %  Tr_Clust_data = Tr_Clust_INS;
                %  TrialType = TrialType_INS;
            elseif strcmp(param.predic{pr}(1),'P') %- pav
                SPK_data = SPK_PAV;
                BEHAV = BEHAV_PAV;
                %   Tr_Clust_data = Tr_Clust_PAV;
                %   TrialType = TrialType_PAV;
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
            bins(~ismember(bins,param.bins4decoding)) = [];

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
                dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' 'rew' 'I_chosenproba' 'I_chosenjuice' }),:)';
                dumm( dumm(:,2)==-1,2)=0;
                if strcmp(param.predic{pr},'I_proba_rew')
                dumm = [dumm(:,1) dumm(:,2)+(10*dumm(:,3))];
                else
                dumm = [dumm(:,1) dumm(:,2)+(10*dumm(:,4))+(10*dumm(:,3))];
                end
            end

            factor = dumm(:,[2 1]);

            %- trials to consider for instrum
            if strcmp(param.predic{pr}(1),'I')
                when = @(arg,cd) BEHAV(strcmp(TrialType_header,arg),:)==cd ;
                diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
                %    diff_juice_all = repmat(diff_juice',length(unique(unit_nb)),1);
            end

            %  predic_all = repmat(dumm(:,2),length(unique(unit_nb)),1);
            %  factor = [dumm(:,2) , unit_nb];

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
            elseif strcmp(param.predic{pr},'I_proba_juice_rew')
                  keep = [310 311 320 321 510 511 520 521 710 711 720 721 910 911 920 921] ;
                remove = ~ismember(dumm(:,2),keep) | ~diff_juice';
            elseif strcmp(param.predic{pr},'I_proba_rew')
                  keep = [300 301 500 501 700 701 900 901 ] ;
                 
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

            conditions_ALL = unique(factor(:,1));
            conditions_rew = unique(mod(factor(:,1),10));
            conditions_pb = unique(round(mod(factor(:,1)/100,100)));
            for n = 1 : length(unit_nb)
                data_sub = data( ismember(factor(:,2),unit_nb(n)),:);
                factor_sub = factor(ismember(factor(:,2),unit_nb(n)),1);

                data_norm = (data_sub-mean(data_sub(:)))/std(data_sub(:));
                data_centered = data_norm-mean(data_norm);

                for f = 1 : length(conditions_ALL)
                    ALLdata{s}(n,f,:) = mean(data_centered(factor_sub==conditions_ALL(f),:));
                end          
                for pp = 1 : 100
                    factor_sub_perm = factor_sub(randperm(length(factor_sub)));
                    for f = 1 : length(conditions_ALL)
                        PERMdata{s,pp}(n,f,:) = mean(data_centered(factor_sub_perm==conditions_ALL(f),:));
                    end      
                end
                for f = 1 : length(conditions_pb)
                    fact = round(mod(factor_sub(:,1)/100,100));
                    PBdata{s}(n,f,:) = mean(data_centered(fact==conditions_pb(f),:));
                end                
                for f = 1 : length(conditions_rew)
                    fact = mod(factor_sub(:,1),10);
                    REWdata{s}(n,f,:) = mean(data_centered(fact==conditions_rew(f),:));
                end
            end
        end
    end

    save('C:\Users\Fred\Documents\POTTc_PCA_subspaces_perm.mat','time','bins','area2test','REWdata','PBdata',"ALLdata",'PERMdata','-v7.3')
end

%% subspace PCA for chosenProba and Reward + Figure 9A

clear
load('C:\Users\Fred\Documents\POTTc_PCA_subspaces_perm.mat')

clear distperm_pb distperm_rew
for ar= 1 : length(area2test)
    for pe = 1 : 100
        neuron_selection =  randperm(size(ALLdata{ar},1),100);

        ALL = [squeeze(ALLdata{ar}(neuron_selection,1,:)) , ...
            squeeze(ALLdata{ar}(neuron_selection,2,:)) , ...
            squeeze(ALLdata{ar}(neuron_selection,3,:)) , ...
            squeeze(ALLdata{ar}(neuron_selection,4,:)) , ...
            squeeze(ALLdata{ar}(neuron_selection,5,:)) , ...
            squeeze(ALLdata{ar}(neuron_selection,6,:)) , ...
            squeeze(ALLdata{ar}(neuron_selection,7,:)) , ...
            squeeze(ALLdata{ar}(neuron_selection,8,:)) ]';

        rejected = find(sum(isnan(ALL))>0);
        %disp(['Removing ' num2str(length(rejected)) ' neurons (FR < ' num2str(thr) ' Hz)'])
        ALL(:,sum(isnan(ALL))>0) = [];
        ALL(:,sum(~isfinite(ALL))>0) = [];

        %- Perform PCA on average FR for each neuron and each proba, at the time of the stim (300 to 460 ms = 26 to 34
        clear PCpb PCrew
        % X = bsxfun(@minus,X,mean(X)); %- no need, it will be done during the PCA

        PCrew.X = nanmean(REWdata{ar}(neuron_selection,:, bins == 6 & (time>100 & time<700)),3)';
        PCpb.X = nanmean(PBdata{ar}(neuron_selection,:, bins == 2 & (time>100 & time<700)),3)';
        % PC.X = mean( all_data_pav(:,:,26:34),3)';
        PCrew.nComp=1;
        PCpb.nComp=1;
        [PCrew.eigenvectors,PCrew.score,PCrew.eigenvalues,~,PCrew.explained,PCrew.mu] = pca(PCrew.X,'NumComponents',PCrew.nComp);
        [PCpb.eigenvectors,PCpb.score,PCpb.eigenvalues,~,PCpb.explained,PCpb.mu] = pca(PCpb.X,'NumComponents',PCpb.nComp);

        PCpb.Xnew = ALL; 
        PCrew.Xnew = ALL; 

        PCrew.Xnew_proj_nonnormalized = (PCrew.Xnew * PCrew.eigenvectors(:,1:PCrew.nComp)) ;
        PCpb.Xnew_proj_nonnormalized = (PCpb.Xnew * PCpb.eigenvectors(:,1:PCpb.nComp)) ;


        pp = randperm(size(PERMdata,2),1); % select 1 permutation

        PERM = [squeeze(PERMdata{ar,pp}(neuron_selection,1,:)) , ...
            squeeze(PERMdata{ar,pp}(neuron_selection,2,:)) , ...
            squeeze(PERMdata{ar,pp}(neuron_selection,3,:)) , ...
            squeeze(PERMdata{ar,pp}(neuron_selection,4,:)) , ...
            squeeze(PERMdata{ar,pp}(neuron_selection,5,:)) , ...
            squeeze(PERMdata{ar,pp}(neuron_selection,6,:)) , ...
            squeeze(PERMdata{ar,pp}(neuron_selection,7,:)) , ...
            squeeze(PERMdata{ar,pp}(neuron_selection,8,:)) ]';

        rejected = find(sum(isnan(PERM))>0);
        %disp(['Removing ' num2str(length(rejected)) ' neurons (FR < ' num2str(thr) ' Hz)'])
        PERM(:,sum(isnan(PERM))>0) = [];
        PERM(:,sum(~isfinite(PERM))>0) = [];

        PCrew.Xperm_proj_nonnormalized = (PERM * PCrew.eigenvectors(:,1:PCrew.nComp)) ;
        PCpb.Xperm_proj_nonnormalized = (PERM * PCpb.eigenvectors(:,1:PCpb.nComp)) ;


        % to reconstruct the time x condition
        findtime = [1:length(time):size(PCrew.Xnew_proj_nonnormalized,1) size(PCrew.Xnew_proj_nonnormalized,1)+1];

        clear score_pb score_rew score_pb_perm score_rew_perm
        for ii = 1 : length(findtime)-1
            score_pb(:,ii) = PCpb.Xnew_proj_nonnormalized(findtime(ii):findtime(ii+1)-1,:);
            score_rew(:,ii) = PCrew.Xnew_proj_nonnormalized(findtime(ii):findtime(ii+1)-1,:);
            score_pb_perm(:,ii) = PCpb.Xperm_proj_nonnormalized(findtime(ii):findtime(ii+1)-1,:);
            score_rew_perm(:,ii) = PCrew.Xperm_proj_nonnormalized(findtime(ii):findtime(ii+1)-1,:);
        end

        %- extract distance
        for t = 1 : length(time)
            dumm_pb = score_pb(t,:);
            dumm_rew = score_rew(t,:);
            dummperm_pb = score_pb_perm(t,:);
            dummperm_rew = score_rew_perm(t,:);
            xx = 0;
            for cd1 = 1 : length(dumm_pb)
                for cd2 = 1 : length(dumm_pb)
                    if cd1 ~= cd2 & cd1>cd2
                        xx = xx + 1;
                        dist_pb{ar}(pe,t,xx) = sqrt(sum( (dumm_pb(cd1) - dumm_pb(cd2))  .^ 2));
                        distperm_pb{ar}(pe,t,xx) = sqrt(sum( (dummperm_pb(cd1) - dummperm_pb(cd2))  .^ 2));
                        dist_rew{ar}(pe,t,xx) = sqrt(sum( (dumm_rew(cd1) - dumm_rew(cd2))  .^ 2));
                        distperm_rew{ar}(pe,t,xx) = sqrt(sum( (dummperm_rew(cd1) - dummperm_rew(cd2))  .^ 2));
                    end
                end
            end
        end

        PCpb_all{ar}.Xnew_proj_nonnormalized(:,pp) = PCpb.Xnew_proj_nonnormalized;
        PCrew_all{ar}.Xnew_proj_nonnormalized(:,pp) = PCrew.Xnew_proj_nonnormalized;
    end

end

figure;
for ar = 1 : length(area2test)
% to reconstruct the time x condition
findtime = [1:length(time):size(PCrew.Xnew_proj_nonnormalized,1) size(PCrew.Xnew_proj_nonnormalized,1)+1];

%- for plotting over time
subsp = 20;
step = 20; % step in time for plotting purposes only (write down value every X time bin)
timesel = [2 4 5 6];
timesel_sub = [-.4 .98 ; 0 .98 ; -.2 .38 ; -.2 .8];

colors = cbrewer('div', 'RdBu', 8);colors=colors([1 8 2 7 3 6 4 5],:);

gaps = [1 find(diff(time)~=subsp)+1 length(time)];
for t = 1 : length(timesel)
    time_considered = find(bins == timesel(t) & time>=1000*timesel_sub(t,1) & time<=1000*timesel_sub(t,2));
    plot_me{t} = time_considered;
end

subplot(3,length(area2test),ar)
timestart = 1;
xlab = [];
col = [1 3 5 7];
% for i = 1 : length(PCpb.score)
%     plot(10,PCpb.score(i),'.','MarkerSize',20,'Color',colors(col(i),:));hold on
% end
for t = 1 : length(timesel)
    timeend = timestart+length(plot_me{t})-1;
    for ii = 1 : length(findtime)-1
        score_curr = nanmean(PCpb_all{ar}.Xnew_proj_nonnormalized(findtime(ii):findtime(ii+1)-1,:),2);
      %  score_curr = PCpb.Xperm_proj_nonnormalized(findtime(ii):findtime(ii+1)-1,:);
        plot(timestart:timeend,score_curr(plot_me{t},1),'Color',colors(ii,:),'LineWidth',2); hold on
    end
    line([timeend timeend],[round(min(PCpb.Xnew_proj_nonnormalized)-1) round(max(PCpb.Xnew_proj_nonnormalized)+1)],'Color',[.6 .6 .6])
    timestart = timeend;
    xlab = [xlab time(plot_me{t})];
end
set(gca,'Xtick',1:step:length(xlab),'XtickLabel',xlab(1:step:end)/1000,'XtickLabelRotation',30,'FontSize',16)
title(area2test{ar})

subplot(3,length(area2test),length(area2test)+ar)
timestart = 1;
xlab = [];
for t = 1 : length(timesel)
    timeend = timestart+length(plot_me{t})-1;
    for ii = 1 : length(findtime)-1
        score_curr = nanmean(PCrew_all{ar}.Xnew_proj_nonnormalized(findtime(ii):findtime(ii+1)-1,:),2);
 %       score_curr = PCrew.Xperm_proj_nonnormalized(findtime(ii):findtime(ii+1)-1,:);
        plot(timestart:timeend,score_curr(plot_me{t},1),'Color',colors(ii,:),'LineWidth',2); hold on
    end
    line([timeend timeend],[round(min(PCrew.Xnew_proj_nonnormalized)-1) round(max(PCrew.Xnew_proj_nonnormalized)+1)],'Color',[.6 .6 .6])
    timestart = timeend;
    xlab = [xlab time(plot_me{t})];
end
set(gca,'Xtick',1:step:length(xlab),'XtickLabel',xlab(1:step:end)/1000,'XtickLabelRotation',30,'FontSize',16)
end

%% Figure 9B

conditions_ALL = [300 301 500 501 700 701 900 901];
%- extract only comparisons between different proba, not same proba but different rew
all_comparisons =[];
xx = 0;
for cd1 = 1 : length(conditions_ALL)
    for cd2 = 1 : length(conditions_ALL)
        if cd1 ~= cd2 & cd1>cd2
            xx = xx + 1;
            all_comparisons(xx,:) =  [conditions_ALL(cd1) , conditions_ALL(cd2)];
        end
    end
end
only_pbdiff = mod(all_comparisons(:,1)-all_comparisons(:,2),10)==0;

%- for plotting over time
subsp = 20;
step = 20; % step in time for plotting purposes only (write down value every X time bin)
timesel = [2 4 5 6];
timesel_sub = [-.4 .98 ; 0 .98 ; -.2 .38 ; -.2 .8];

gaps = [1 find(diff(time)~=subsp)+1 length(time)];
for t = 1 : length(timesel)
    time_considered = find(bins == timesel(t) & time>=1000*timesel_sub(t,1) & time<=1000*timesel_sub(t,2));
    plot_me{t} = time_considered;
end

colorsPB = cbrewer('seq', 'Greens', 6);
colorsREW = cbrewer('seq', 'Greys', 6);

figure;
for ar = 1 : length(area2test)

    dist_m1 = mean(dist_pb{ar},3);
    distperm_m1 = mean(distperm_pb{ar},3);
    dist_m1b = mean(dist_pb{ar}(:,:,only_pbdiff),3);
    distperm_m1b = mean(distperm_pb{ar}(:,:,only_pbdiff),3);
    dist_m2 = mean(dist_rew{ar},3);
    distperm_m2 = mean(distperm_rew{ar},3);

    %- different way of checking stats
    clear pval_pb_test pval_rew_test
    pval_pb_test_sig = false(size(dist_m1,2),size(dist_m1,1));
    pval_rew_test_sig = false(size(dist_m2,2),size(dist_m2,1));
    for tt = 1 : size(dist_m1,2)
        for nn = 1 : size(dist_m1,1)
            pval_pb_test(tt,nn) = sum(distperm_m1(:,tt)>dist_m1(nn,tt)) /length(distperm_m1(:,tt));
            pval_rew_test(tt,nn) = sum(distperm_m2(:,tt)>dist_m2(nn,tt)) /length(distperm_m2(:,tt));
        end
            [~,idxs_pb]=findenough(pval_pb_test(tt,:),0.01,4,'<');
            [~,idxs_rew]=findenough(pval_rew_test(tt,:),0.01,4,'<');
        if ~isempty(idxs_pb)
            pval_pb_test_sig(tt,idxs_pb) = true;
        end
        if ~isempty(idxs_rew)
            pval_rew_test_sig(tt,idxs_rew) = true;
        end
    end
    
    showsig_pb = NaN(size(pval_pb_test));
    showsig_pb((100*mean(pval_pb_test_sig,2))>50) = true;
    showsig_rew = NaN(size(pval_rew_test));
    showsig_rew((100*mean(pval_rew_test_sig,2))>50) = true;

    subplot(3,3,ar)

    timestart = 1;
    for t = 1 : length(timesel)
        timeend = timestart+length(plot_me{t})-1;

        plot(timestart:timeend,mean(dist_m1b(:,plot_me{t})),'--','Color',colorsPB(5,:),'LineWidth',2); hold on
       
        plot(timestart:timeend,mean(dist_m1(:,plot_me{t})),'Color',colorsPB(5,:),'LineWidth',2); hold on
        ciplot(mean(dist_m1(:,plot_me{t}))-std(dist_m1(:,plot_me{t})),mean(dist_m1(:,plot_me{t}))+std(dist_m1(:,plot_me{t})),timestart:timeend,colorsPB(4,:),1);
        plot(timestart:timeend,mean(distperm_m1(:,plot_me{t})),'Color',colorsPB(5,:),'LineWidth',1); hold on
        plot(timestart:timeend,5.6*showsig_pb(plot_me{t}),'.','Color',colorsPB(5,:))
        %  plot(timestart:timeend,mean(dist_m1(:,plot_me{t}))+std(dist_m1(:,plot_me{t})),'Color',colorsPB(3,:),'LineWidth',2); hold on
        %  plot(timestart:timeend,mean(dist_m1(:,plot_me{t}))-std(dist_m1(:,plot_me{t})),'Color',colorsPB(3,:),'LineWidth',2); hold on
        plot(timestart:timeend,mean(dist_m2(:,plot_me{t})),'Color',colorsREW(5,:),'LineWidth',2); hold on
        plot(timestart:timeend,mean(distperm_m2(:,plot_me{t})),'Color',colorsREW(5,:),'LineWidth',1); hold on
        ciplot(mean(dist_m2(:,plot_me{t}))-std(dist_m2(:,plot_me{t})),mean(dist_m2(:,plot_me{t}))+std(dist_m2(:,plot_me{t})),timestart:timeend,colorsREW(4,:),1);
        plot(timestart:timeend,5.7*showsig_rew(plot_me{t}),'.','Color',colorsREW(5,:))
        %   plot(timestart:timeend,mean(dist_m2(:,plot_me{t}))+std(dist_m2(:,plot_me{t})),'Color',colorsREW(3,:),'LineWidth',2); hold on
        %   plot(timestart:timeend,mean(dist_m2(:,plot_me{t}))-std(dist_m2(:,plot_me{t})),'Color',colorsREW(3,:),'LineWidth',2); hold on
        timestart = timeend;
        xlab = [xlab time(plot_me{t})];
        ylim([0 6])
    end
    set(gca,'Xtick',1:step:length(xlab),'XtickLabel',xlab(1:step:end)/1000,'XtickLabelRotation',30,'FontSize',16)
    title(area2test{ar})
end
