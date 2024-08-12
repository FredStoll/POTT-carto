%% LDA on pseudo-population for Stimulus and Reward period
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.04

%% STIMULUS

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
list = dir([path2go 'POOL_2*.mat']);
% list = dir([path2go 'POOL_10*.mat']);
% list = dir([path2go 'POOL_18*.mat']);
params = load([path2go 'M021519a_SPKpool.mat']);

param.predic = {'I_chosenjuice' 'I_chosenproba' 'I_chosenside'}; %- select param you want to test..
param.minNeurons = [25 50 75 100 125 150 175 200 300 400 500 600 700 800 900 1000]; %- min number of neuron to run
param.nComp = 20; %- nb of component to keep for PCA
param.Repetition = 200; %- number of times to run the decoder (100 or 200 ideally)
param.minTrCond = 10; %- minimum number of trials per condition to run
param.overwrite = false;
param.thresh = false;
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
if exist([path2go 'res_LDApop_kfold_final_' list(1).name(un(1)+1:un(2)-1) '.mat'])==2
    prev = load([path2go 'res_LDApop_kfold_final_' list(1).name(un(1)+1:un(2)-1) '.mat']);
    done = [];
    for pr = 1 : length(param.predic)
        if isfield(prev.res_LDApop,param.predic{pr})
            done(pr) = true;
        end
    end
    if ~param.overwrite
        param.predic(find(done==true))=[];
    end
    res_LDApop = prev.res_LDApop;
end
if param.overwrite
    clear res_LDApop
end

if ~isempty(param.predic) %- if there is some predictors not done yet
    
    for pr = 1 : length(param.predic)
        %nSess = zeros(length(area2test),1);
        for s = 1 : length(list)
            clearvars -except path2go list param params areas area2test util s u pr res_LDApop nSess un
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
                for n = 1 : length(neurons_rank)
                    if mean(mean(SPK_data_cut(Tr_Clust_data(:,2)==n,:)))<param.thr  %- put NaNs when FR is below thr
                        reject = [reject , n];
                    end
                end
                
                %- remove them from main matrices
                SPK_data_cut( ismember(Tr_Clust_data(:,2),reject) ,:) = [];
                Tr_Clust_data( ismember(Tr_Clust_data(:,2),reject) ,:) = [];
                area_histology(reject)=[];
            end
            
            %- create data and factor matrix for the decoding
            data = SPK_data_cut;
            unit_nb = unique(BEHAV(1,:));
            
            TrialType_header = ['Unit_Nb' 'Session_Nb' params.TrialType_header];
            
            dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' param.predic{pr} }),:)';
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
                        
                        %- different selection of neurons
                        %   unit_nb = unique(BEHAV(1,:));
                        %   unit_nb( sum(nTr>param.minTrCond,2)~= length(keep) ) = [];
                        
                        if ~isempty(param.minNeurons(u))
                            if param.minNeurons(u)<=length(unit_nb)
                                rd = randperm(length(unit_nb));
                                unit_nb_sub = unit_nb(rd);
                                unit_nb_sub = unit_nb_sub(1:param.minNeurons(u));
                            end
                        end
                        
                        data_sub =  data( ismember(factor(:,2),unit_nb_sub),:);
                        factor_sub =  factor( ismember(factor(:,2),unit_nb_sub),:);
                        
                        [XX,YY,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,1,'perm',1,'minTr',[],'pop',false); % ,'window',[-.5 1.5]
                        [perf(:,p),~] = pca_lda_kfold(XX,YY,param);
                        
                        %- permute trial factors
                        if param.tr_perm
                            loopi = unique(factor_sub(:,2));
                            for up = 1 : length(loopi)
                                dumm_condi = factor_sub(factor_sub(:,2)==loopi(up),1);
                                
                                factor_sub(factor_sub(:,2)==loopi(up),1) = dumm_condi(randperm(length(dumm_condi)));
                            end
                            
                            [XXperm,YYperm,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,1,'perm',1,'minTr',[],'pop',false); % ,'window',[-.5 1.5]
                            [perf_perm(:,p),~] = pca_lda_kfold(XXperm,YYperm,param);
                        end
                    end
                    
                    %- extract param of interest to save
                    %  nSess(ar,1) = nSess(ar,1) + 1 ;
                    temp.lda_sess = name;
                    temp.perf = perf;
                    temp.nbUnit = param.minNeurons(u);
                    if param.tr_perm
                        temp.perf_perm = perf_perm;
                    end
                    unde = find(name=='_',1,'last');
                    area_name = ['area_' name(unde+1:end-4)];
                    eval(['res_LDApop.' param.predic{pr} '.' area_name '(u) = temp;'])
                    
                end
            end
        end
        eval(['res_LDApop.' param.predic{pr} '.param = param;'])
    end
    
    save([path2go 'res_LDApop_kfold_final_' list(1).name(un(1)+1:un(2)-1) '.mat'],'res_LDApop','list','area2test','path2go')
end


%% REWARD

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
list = dir([path2go 'POOL_2*.mat']);
% list = dir([path2go 'POOL_10*.mat']);
% list = dir([path2go 'POOL_18*.mat']);
params = load([path2go 'M021519a_SPKpool.mat']);

param.predic = {'I_rew' 'I_chosenjuicerew'}; %- select param you want to test..
param.minNeurons = [25 50 75 100 125 150 175 200 300 400 500 600 700 800 900 1000]; %- min number of neuron to run
param.nComp = 20; %- nb of component to keep for PCA
param.Repetition = 200; %- number of times to run the decoder (100 or 200 ideally)
param.minTrCond = 10; %- minimum number of trials per condition to run
param.overwrite = false;
param.thresh = false;
param.thr = 1; % in Hz, used if param.thresh is 'y'
param.rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
param.bins4decoding = [6] ; %- 2 // 6 for rew %- perform decoding on subset on bins (stim period here)
param.timebin = [100 700] ;% 200 800 // 0 600 rew %- perform decoding on subset on bins (stim period here)
param.normalize_fr = false;
param.tr_perm = true;

disp(['Computing LDA on ' num2str(length(list)) ' sessions'])
areas = utils_POTT_areas;

area2test = {'12r' '12m' '12o' '12l' 'a11ml' '13l' '13m' 'LAI'};

util.minmaxnorm = @(data) (data-min(data))/(max(data)-min(data));
util.meannorm  = @(data) (data-mean(data))/(max(data)-min(data));
util.minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

un = find(list(1).name=='_');

%- check if matrix already exists for each predictor
if exist([path2go 'res_LDApop_rew_kfold_v2_' list(1).name(un(1)+1:un(2)-1) '.mat'])==2
    prev = load([path2go 'res_LDApop_rew_kfold_v2_' list(1).name(un(1)+1:un(2)-1) '.mat']);
    done = [];
    for pr = 1 : length(param.predic)
        if isfield(prev.res_LDApop,param.predic{pr})
            done(pr) = true;
        end
    end
    if ~param.overwrite
        param.predic(find(done==true))=[];
    end
    res_LDApop = prev.res_LDApop;
end
if param.overwrite
    clear res_LDApop
end

if ~isempty(param.predic) %- if there is some predictors not done yet
    
    for pr = 1 : length(param.predic)
        %nSess = zeros(length(area2test),1);
        for s = 1 : length(list)
            clearvars -except path2go list param params areas area2test util s u pr res_LDApop nSess un
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
                for n = 1 : length(neurons_rank)
                    if mean(mean(SPK_data_cut(Tr_Clust_data(:,2)==n,:)))<param.thr  %- put NaNs when FR is below thr
                        reject = [reject , n];
                    end
                end
                
                %- remove them from main matrices
                SPK_data_cut( ismember(Tr_Clust_data(:,2),reject) ,:) = [];
                Tr_Clust_data( ismember(Tr_Clust_data(:,2),reject) ,:) = [];
                area_histology(reject)=[];
            end
            
            %- create data and factor matrix for the decoding
            data = SPK_data_cut;
            unit_nb = unique(BEHAV(1,:));
            
            TrialType_header = ['Unit_Nb' 'Session_Nb' params.TrialType_header];
            
            if strcmp(param.predic{pr},'P_rew') | strcmp(param.predic{pr},'I_rew') 
               dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' 'rew' }),:)';
            elseif strcmp(param.predic{pr},'P_juicerew')   
               dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' 'P_juice' }),:)';
            elseif  strcmp(param.predic{pr},'I_chosenjuicerew')
               dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' 'I_chosenjuice' }),:)';
            else
                dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' param.predic{pr} }),:)';
            end
            factor = dumm(:,[2 1]);
            
            %- trials to consider for instrum
            if strcmp(param.predic{pr}(1),'I')
                when = @(arg,cd) BEHAV(strcmp(TrialType_header,arg),:)==cd ;
                diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
            end

            if strcmp(param.predic{pr}(end-2:end),'rew')
                when = @(arg,cd) BEHAV(strcmp(TrialType_header,arg),:)==cd ;
                rewarded = when('rew',1) ; %- take only rewarded trials
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

            elseif strcmp(param.predic{pr},'P_juicerew')   
                keep = [1 2];
                remove = ~rewarded';  % was not removing same juice trials before
            elseif  strcmp(param.predic{pr},'I_chosenjuicerew')
                keep = [1 2];
                remove = ~diff_juice' & ~rewarded';  % was not removing same juice trials before
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
                        
                        %- different selection of neurons
                        %   unit_nb = unique(BEHAV(1,:));
                        %   unit_nb( sum(nTr>param.minTrCond,2)~= length(keep) ) = [];
                        
                        if ~isempty(param.minNeurons(u))
                            if param.minNeurons(u)<=length(unit_nb)
                                rd = randperm(length(unit_nb));
                                unit_nb_sub = unit_nb(rd);
                                unit_nb_sub = unit_nb_sub(1:param.minNeurons(u));
                            end
                        end
                        
                        data_sub =  data( ismember(factor(:,2),unit_nb_sub),:);
                        factor_sub =  factor( ismember(factor(:,2),unit_nb_sub),:);
                        
                        [XX,YY,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,1,'perm',1,'minTr',[],'pop',false); % ,'window',[-.5 1.5]
                        [perf(:,p),~] = pca_lda_kfold(XX,YY,param);
                        
                        %- permute trial factors
                        if param.tr_perm
                            loopi = unique(factor_sub(:,2));
                            for up = 1 : length(loopi)
                                dumm_condi = factor_sub(factor_sub(:,2)==loopi(up),1);
                                
                                factor_sub(factor_sub(:,2)==loopi(up),1) = dumm_condi(randperm(length(dumm_condi)));
                            end
                            
                            [XXperm,YYperm,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,1,'perm',1,'minTr',[],'pop',false); % ,'window',[-.5 1.5]
                            [perf_perm(:,p),~] = pca_lda_kfold(XXperm,YYperm,param);
                        end
                    end
                    
                    %- extract param of interest to save
                    %  nSess(ar,1) = nSess(ar,1) + 1 ;
                    temp.lda_sess = name;
                    temp.perf = perf;
                    temp.nbUnit = param.minNeurons(u);
                    if param.tr_perm
                        temp.perf_perm = perf_perm;
                    end
                    unde = find(name=='_',1,'last');
                    area_name = ['area_' name(unde+1:end-4)];
                    eval(['res_LDApop.' param.predic{pr} '.' area_name '(u) = temp;'])
                end
            end
        end
        eval(['res_LDApop.' param.predic{pr} '.param = param;'])
    end
    
    save([path2go 'res_LDApop_rew_kfold_final_' list(1).name(un(1)+1:un(2)-1) '.mat'],'res_LDApop','list','area2test','path2go')
end


%% REVISION = ChosenFlavor in unrewarded trials


clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'; %- path where SPKpool files are!
% list = dir([path2go 'POOL_2*.mat']);
% list = dir([path2go 'POOL_10*.mat']);
 list = dir([path2go 'POOL_18*.mat']);
params = load([path2go 'M021519a_SPKpool.mat']);

param.predic = {'I_chosenjuicerew'}; %- select param you want to test..
param.minNeurons = [25 50 75 100 125 150 175 200 300 400 500 600 700 800 900 1000]; %- min number of neuron to run
param.nComp = 20; %- nb of component to keep for PCA
param.Repetition = 200; %- number of times to run the decoder (100 or 200 ideally)
param.minTrCond = 10; %- minimum number of trials per condition to run
param.overwrite = false;
param.thresh = false;
param.thr = 1; % in Hz, used if param.thresh is 'y'
param.rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
param.bins4decoding = [6] ; %- 2 // 6 for rew %- perform decoding on subset on bins (stim period here)
param.timebin = [100 700] ;% 200 800 // 0 600 rew %- perform decoding on subset on bins (stim period here)
param.normalize_fr = false;
param.tr_perm = true;

disp(['Computing LDA on ' num2str(length(list)) ' sessions'])
areas = utils_POTT_areas;

area2test = {'12r' '12m' '12o' '12l' 'a11ml' '13l' '13m' 'LAI'};

util.minmaxnorm = @(data) (data-min(data))/(max(data)-min(data));
util.meannorm  = @(data) (data-mean(data))/(max(data)-min(data));
util.minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

un = find(list(1).name=='_');

%- check if matrix already exists for each predictor
if exist([path2go 'res_LDApop_rew_kfold_v2_' list(1).name(un(1)+1:un(2)-1) '.mat'])==2
    prev = load([path2go 'res_LDApop_rew_kfold_v2_' list(1).name(un(1)+1:un(2)-1) '.mat']);
    done = [];
    for pr = 1 : length(param.predic)
        if isfield(prev.res_LDApop,param.predic{pr})
            done(pr) = true;
        end
    end
    if ~param.overwrite
        param.predic(find(done==true))=[];
    end
    res_LDApop = prev.res_LDApop;
end
if param.overwrite
    clear res_LDApop
end

if ~isempty(param.predic) %- if there is some predictors not done yet
    
    for pr = 1 : length(param.predic)
        %nSess = zeros(length(area2test),1);
        for s = 1 : length(list)
            clearvars -except path2go list param params areas area2test util s u pr res_LDApop nSess un
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
                for n = 1 : length(neurons_rank)
                    if mean(mean(SPK_data_cut(Tr_Clust_data(:,2)==n,:)))<param.thr  %- put NaNs when FR is below thr
                        reject = [reject , n];
                    end
                end
                
                %- remove them from main matrices
                SPK_data_cut( ismember(Tr_Clust_data(:,2),reject) ,:) = [];
                Tr_Clust_data( ismember(Tr_Clust_data(:,2),reject) ,:) = [];
                area_histology(reject)=[];
            end
            
            %- create data and factor matrix for the decoding
            data = SPK_data_cut;
            unit_nb = unique(BEHAV(1,:));
            
            TrialType_header = ['Unit_Nb' 'Session_Nb' params.TrialType_header];
            
            if strcmp(param.predic{pr},'P_rew') | strcmp(param.predic{pr},'I_rew') 
               dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' 'rew' }),:)';
            elseif strcmp(param.predic{pr},'P_juicerew')   
               dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' 'P_juice' }),:)';
            elseif  strcmp(param.predic{pr},'I_chosenjuicerew')
               dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' 'I_chosenjuice' }),:)';
            else
                dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' param.predic{pr} }),:)';
            end
            factor = dumm(:,[2 1]);
            
            %- trials to consider for instrum
            if strcmp(param.predic{pr}(1),'I')
                when = @(arg,cd) BEHAV(strcmp(TrialType_header,arg),:)==cd ;
                diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
            end

            if strcmp(param.predic{pr}(end-2:end),'rew')
                when = @(arg,cd) BEHAV(strcmp(TrialType_header,arg),:)==cd ;
                rewarded = when('rew',-1) ; %- take only rewarded trials %%%%%%%%%% THIS IS WHAT CHANGED!!!!
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

            elseif strcmp(param.predic{pr},'P_juicerew')   
                keep = [1 2];
                remove = ~rewarded';  % was not removing same juice trials before
            elseif  strcmp(param.predic{pr},'I_chosenjuicerew')
                keep = [1 2];
                remove = ~diff_juice' & ~rewarded';  % was not removing same juice trials before
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
                        
                        %- different selection of neurons
                        %   unit_nb = unique(BEHAV(1,:));
                        %   unit_nb( sum(nTr>param.minTrCond,2)~= length(keep) ) = [];
                        
                        if ~isempty(param.minNeurons(u))
                            if param.minNeurons(u)<=length(unit_nb)
                                rd = randperm(length(unit_nb));
                                unit_nb_sub = unit_nb(rd);
                                unit_nb_sub = unit_nb_sub(1:param.minNeurons(u));
                            end
                        end
                        
                        data_sub =  data( ismember(factor(:,2),unit_nb_sub),:);
                        factor_sub =  factor( ismember(factor(:,2),unit_nb_sub),:);
                        
                        [XX,YY,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,1,'perm',1,'minTr',[],'pop',false); % ,'window',[-.5 1.5]
                        [perf(:,p),~] = pca_lda_kfold(XX,YY,param);
                        
                        %- permute trial factors
                        if param.tr_perm
                            loopi = unique(factor_sub(:,2));
                            for up = 1 : length(loopi)
                                dumm_condi = factor_sub(factor_sub(:,2)==loopi(up),1);
                                
                                factor_sub(factor_sub(:,2)==loopi(up),1) = dumm_condi(randperm(length(dumm_condi)));
                            end
                            
                            [XXperm,YYperm,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,1,'perm',1,'minTr',[],'pop',false); % ,'window',[-.5 1.5]
                            [perf_perm(:,p),~] = pca_lda_kfold(XXperm,YYperm,param);
                        end
                    end
                    
                    %- extract param of interest to save
                    %  nSess(ar,1) = nSess(ar,1) + 1 ;
                    temp.lda_sess = name;
                    temp.perf = perf;
                    temp.nbUnit = param.minNeurons(u);
                    if param.tr_perm
                        temp.perf_perm = perf_perm;
                    end
                    unde = find(name=='_',1,'last');
                    area_name = ['area_' name(unde+1:end-4)];
                    eval(['res_LDApop.' param.predic{pr} '.' area_name '(u) = temp;'])
                end
            end
        end
        eval(['res_LDApop.' param.predic{pr} '.param = param;'])
    end
    
    save([path2go 'res_LDApop_unrew_kfold_final_' list(1).name(un(1)+1:un(2)-1) '.mat'],'res_LDApop','list','area2test','path2go')
end

