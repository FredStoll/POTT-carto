%% ANOVA for POTT dataset

%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.03

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/neurons/subset-carto/'
list = dir([path2go '*a_SPKpool.mat']);

disp(['Computing ANOVANs on ' num2str(length(list)) ' sessions'])

rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
%- times_evts = {'FixFP_onset' 'Stim_onset' 'Resp_onset' 'FixResp_onset' 'FB_onset' 'Rew_onset' 'FB_offset'};
bins2consider=[1 2 3 4 5 6]; %- perform analyses on subset of bins (relates too times_evts)

%- normalization
normalize_fr = true;
minmaxnorm = @(data) (data-min(data))/(max(data)-min(data));
meannorm  = @(data) (data-mean(data))/(max(data)-min(data));
norm_pb = @(data) -1+((data-0.1)*2)/(0.9-0.1) ;
minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

predic={'I_chosenproba' ; 'I_unchosenproba' ; 'I_chosenjuice' ; 'I_chosenside'}; %- select any param you want to test.. order matters!
predic_REW={'rew'; 'I_chosenjuicerew' }; %- select any param you want to test.. order matters!

for s = 1 : length(list)
    name = list(s).name;

    load([path2go name]); disp(s)

    %- remove unwanted bins
    bins2remove = (rmv/subsp)-1;
    n_evt = length(times_evts);
    lin=size(neurons_rank,1)*n_evt;
    nTimes = [0 ; (sum(abs(times(:,:)),2)/subsp)];
    col = sum(nTimes);
    bins = NaN(1,col);
    time = NaN(1,col);

    for b = 1 : length(times_evts)
        time(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) =  (times(b,1):subsp:times(b,2)-(subsp/2));
        bins(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) = b;
        bins([sum(nTimes(1:b))+1:sum(nTimes(1:b))+1+bins2remove   sum(nTimes(1:b+1))-bins2remove:sum(nTimes(1:b+1))    ]) = NaN; %- remove 100 ms each side (avoid overlaps between bins and smoothing problems)
    end

    if ~isempty(SPK_INS)

        %- normalize FR
        if normalize_fr %- normalize FR using min-max in considered bins only (all trials except border of events)
            SPK_data_norm = SPK_INS;
            units = unique(Tr_Clust_INS(:,2));
            for u = 1 : length(units)
                temp = SPK_INS(Tr_Clust_INS(:,2)==units(u),~isnan(bins)); %- could normalize only on the bin2use, but more likely to get some 0 FR?? so NaNs!
                SPK_data_norm(Tr_Clust_INS(:,2)==units(u),~isnan(bins)) = reshape(minnorm(temp(:)),size(temp,1),size(temp,2));
            end

            SPK_INS = SPK_data_norm;
        end

        %- cut out the overlapping bins
        SPK_INS_cut = SPK_INS;
        SPK_INS_cut(:,~ismember(bins,bins2consider))=[];
        time(~ismember(bins,bins2consider)) = [];
        bins_considered = bins(ismember(bins,bins2consider));

        %- create factor matrix given predictors selected
        for pa = 1:length(predic)
            factors_INS(:,pa) = TrialType_INS(ismember(TrialType_header,predic(pa)),:)';
        end

        %- trials to considerfor instrum
        when = @(arg,cd) TrialType_INS(strcmp(TrialType_header,arg),:)==cd ;
        diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
        same_juice = (when('I_juiceL',1) & when('I_juiceR',1)) | (when('I_juiceL',2) & when('I_juiceR',2)); %- take only the different juice trials

        % nbTR(s) = sum(diff_juice' &  ismember(factors_INS(:,1),[10:20:90]) )
        tr2consider = find(diff_juice' &  ismember(factors_INS(:,1),[10:20:90]));

        %- normalize predictors
        factors_INS(:,1) = norm_pb(factors_INS(:,1)/100); %- normalize proba
        factors_INS(:,2) = norm_pb(factors_INS(:,2)/100); %- normalize proba
        factors_INS(factors_INS(:,3)==1,3)=-1; %- standardize it (-1 vs 1)
        factors_INS(factors_INS(:,3)==2,3)=1; 
        factors_INS(factors_INS(:,4)==1,4)=-1; %-  standardize it (-1 vs 1)
        factors_INS(factors_INS(:,4)==2,4)=1; 

        %  if length(tr2consider)>=100
        % tr2consider = tr2consider(1:150);

        parfor n = 1 : size(neurons_rank,1)
            disp(['Neuron ' num2str(n) '/' num2str(size(neurons_rank,1))])

            %- INSTRUMENTAL TASK - Perform ANOVANs
            data = SPK_INS_cut(Tr_Clust_INS(:,2)==neurons_rank(n,1),:);

            for j = 1 : size(data,2) %- perform the ANOVANs for every time bins
                [ANOVA_INS(n).Ap(:,j),At,Amdl,~] = anovan(data(tr2consider,j),factors_INS(tr2consider,[1 3 4]),'continuous',[1],'varname',predic([1 3 4]),'model',diag(ones(length(predic([1 3 4])),1)),'display','off');
                for m = 1 : length(predic([1 3 4]))
                    ANOVA_INS(n).PEVs(m,j) = ((At{m+1,2}/At{end,2}))*100;
                    ANOVA_INS(n).Fs(m,j) = At{m+1,6};
                    ANOVA_INS(n).Omega2(m,j) = (At{m+1,2}-At{m+1,3}*At{end-1,5})/(At{end-1,5}+At{end,2});
                    ANOVA_INS(n).Betas(m,j) = Amdl.coeffs(m+1);
                end
                ANOVA_INS(n).R2(j) = 1 - (At{end-1,2}/At{end,2});
            end
            ANOVA_INS(n).VIF=vif(factors_INS(tr2consider,:));

        end

        for n = 1 : size(neurons_rank,1)
            res_anova(s).ins_diff.Ap(n,:,:) = ANOVA_INS(n).Ap;
            res_anova(s).ins_diff.PEVs(n,:,:) = ANOVA_INS(n).PEVs;
            res_anova(s).ins_diff.Fs(n,:,:) = ANOVA_INS(n).Fs;
            res_anova(s).ins_diff.Omega2(n,:,:) = ANOVA_INS(n).Omega2;
            res_anova(s).ins_diff.Betas(n,:,:) = ANOVA_INS(n).Betas;
            res_anova(s).ins_diff.R2(n,:) = ANOVA_INS(n).R2;
            res_anova(s).ins_diff.VIF = ANOVA_INS(1).VIF;
            res_anova(s).ins_diff.predic = predic([1 3 4]);
        end
        res_anova(s).neurons_area = neurons_area ;
        res_anova(s).neurons_info = neurons_info ;
        res_anova(s).neurons_info_header = neurons_info_header ;
        res_anova(s).session = name(1:8);
        %  end

        %- SECOND MODEL
        %- create factor matrix given predictors selected
        for pa = 1:length(predic_REW)
            if strcmp(predic_REW(pa),'I_chosenjuicerew')
                factors_REW(:,pa) = TrialType_INS(ismember(TrialType_header,'I_chosenjuice'),:)';
            else
                factors_REW(:,pa) = TrialType_INS(ismember(TrialType_header,predic_REW(pa)),:)';
            end
        end
        %- trials to considerfor instrum
        when = @(arg,cd) TrialType_INS(strcmp(TrialType_header,arg),:)==cd ;
        diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
        same_juice = (when('I_juiceL',1) & when('I_juiceR',1)) | (when('I_juiceL',2) & when('I_juiceR',2)); %- take only the different juice trials

        % nbTR(s) = sum(diff_juice' &  ismember(factors_INS(:,1),[10:20:90]) )
        tr2consider = find(diff_juice');

        %- normalize predictors
        factors_REW(factors_REW(:,2)==1,2)=-1; %- a way to standardize it (-1 vs 1)
        factors_REW(factors_REW(:,2)==2,2)=1; %- a way to standardize it (-1 vs 1)

        factors_REW_nested = factors_REW;
        factors_REW_nested(factors_REW_nested(:,1)==-1,2) = 0;

        parfor n = 1 : size(neurons_rank,1)
            disp(['Neuron ' num2str(n) '/' num2str(size(neurons_rank,1))])

            %- INSTRUMENTAL TASK - Perform ANOVANs
            data = SPK_INS_cut(Tr_Clust_INS(:,2)==neurons_rank(n,1),:);

            for j = 1 : size(data,2) %- perform the ANOVANs for every time bins
                [ANOVA_REW(n).Ap(:,j),At,Amdl,~] = anovan(data(tr2consider,j),factors_REW_nested(tr2consider,:),'continuous',[],'nested',[0 0 ; 1 0],'varname',predic_REW,'model',[1 0 ; 0 1],'display','off');
                for m = 1 : length(predic_REW)
                    ANOVA_REW(n).PEVs(m,j) = ((At{m+1,2}/At{end,2}))*100;
                    ANOVA_REW(n).Fs(m,j) = At{m+1,6};
                    ANOVA_REW(n).Omega2(m,j) = (At{m+1,2}-At{m+1,3}*At{end-1,5})/(At{end-1,5}+At{end,2});
                    ANOVA_REW(n).Betas(m,j) = Amdl.coeffs(m+1);
                end
                ANOVA_REW(n).R2(j) = 1 - (At{end-1,2}/At{end,2});
            end
            ANOVA_REW(n).VIF=vif(factors_REW(tr2consider,:));
        end

        for n = 1 : size(neurons_rank,1)
            res_anova(s).rew_diff.Ap(n,:,:) = ANOVA_REW(n).Ap;
            res_anova(s).rew_diff.PEVs(n,:,:) = ANOVA_REW(n).PEVs;
            res_anova(s).rew_diff.Fs(n,:,:) = ANOVA_REW(n).Fs;
            res_anova(s).rew_diff.Omega2(n,:,:) = ANOVA_REW(n).Omega2;
            res_anova(s).rew_diff.Betas(n,:,:) = ANOVA_REW(n).Betas;
            res_anova(s).rew_diff.R2(n,:) = ANOVA_REW(n).R2;
            res_anova(s).rew_diff.VIF = ANOVA_REW(1).VIF;
            res_anova(s).rew_diff.predic = predic_REW;
        end
        res_anova(s).neurons_area = neurons_area ;
        res_anova(s).neurons_info = neurons_info ;
        res_anova(s).neurons_info_header = neurons_info_header ;
        res_anova(s).session = name(1:8);

        clearvars -except minnorm norm_pb res_anova path2go thr rmv s list normalize_fr minmaxnorm meannorm bins2consider bins_considered predic_REW predic predic_PAV_FB overwrite thr_len thresh time subsp times_evts

    end
end
save([path2go 'res_ANOVA_full_nested.mat'],'res_anova','list','bins2consider','bins_considered','time','subsp','times_evts','-v7.3')
