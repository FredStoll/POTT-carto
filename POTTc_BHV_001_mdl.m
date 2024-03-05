%% Main script for BEHAV analysis of the POTTc dataset
%- 
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2022.12

clear

%- locate the files
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-Carto/data/behav/';

skip = true; %- re-run the models or just run on saved models! 

if ~skip
    mk = 'Morbier'; % 'Morbier' 'Mimic'

    %- list all the files
    if strcmp(mk(1:2),'Mi')
        list = dir([path2go 'X*a_behav.mat']);
    else
        list = dir([path2go 'M*a_behav.mat']);
    end

    showfig = false;

    %- reorganize list by date
    days = [];
    for ss = 1 : length(list)
        days = [days;datenum(list(ss).name(2:7),'mmddyy') , ss];
    end
    date_order = sortrows(days);
    list = list(date_order(:,2));

    x=0;
    for currsess = 1 : length(list) % for every sessions
        filename = [path2go list(currsess).name];

        [all,par] = behav_model_proba(filename,false,showfig); %- no permutations
        if ~isempty(all)
            x = x + 1;
            ALL(x,1)=all;
            param(x,1) = par;
        end

    end

    save([path2go mk '_behav_bins_carto.mat'],'ALL','param','-v7.3')
end


%% POTT - load models for behavioral analyses and Figures (1B/4A-C/5)

M1 = load([path2go 'Morbier_behav_bins_carto.mat'],'ALL','param')
M2 = load([path2go 'Mimic_behav_bins_carto.mat'],'ALL','param')

%- remove sessions before threshold for Morbier, if needed
subset = {'052318' '010122'} ; %- for anything after last change
clear days 
if ~isempty(subset)
    lim = datenum(subset,'mmddyy');
    for ss = 1 : length(M1.ALL)
        days(ss) = datenum(M1.ALL(ss).name(end-16:end-11),'mmddyy');
    end
    takeme = (days>=lim(1) & days<lim(2)) ;
else
    takeme = true(1,length(M1.ALL));
end
M1.ALL = M1.ALL(takeme);
M1.param = M1.param(takeme);

%- number trials per sessions (all set, not only converging ones)
x=0;
nTr_tot=[];
for m = 1 : 2
    eval(['ALL = M' num2str(m) '.ALL;']);
    for s = 1 : length(ALL)
        x = x + 1;
        nTr_tot(x,:) = [size(ALL(s).T,1) m];
    end
end
[min(nTr_tot(nTr_tot(:,2)==1,1)) median(nTr_tot(nTr_tot(:,2)==1,1)) max(nTr_tot(nTr_tot(:,2)==1,1)) ]
[min(nTr_tot(nTr_tot(:,2)==2,1)) median(nTr_tot(nTr_tot(:,2)==2,1)) max(nTr_tot(nTr_tot(:,2)==2,1)) ]

all_converge = [M1.ALL(:).converge M2.ALL(:).converge ];

%- save a matrix with the names of the considered sessions + converge model
clear behav_sess
x=0;
for m = 1 : 2
    clear conv
    eval(['ALL = M' num2str(m) '.ALL;']);
   
    takeme = [ALL(:).converge]==1;
    sessions = cat(1,ALL(takeme).name);
    for s = 1 : size(sessions,1)
        x = x + 1;
        behav_sess{x,1} = sessions(s,end-17:end-10);
    end
end

%% Fig 1B-C - Behav model R squared + Choice probability
figure;
colors = cbrewer('qual', 'Paired', 10);

subplot(1,8,1)
%- plot adjusted R2
for m = 1 : 2
    clear converge r_adj  keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    for i = 1 : length(ALL)
        r_adj(:,i) = ALL(i).mdl.Rsquared.Adjusted;
        converge(:,i) = ALL(i).converge;
    end
    keep = converge==1 ;
    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
    else
        col = colors(9:10,:);
        mm = .2;
    end
    X= r_adj(1,keep);
    yl=1+mm;
    wdth = .35;
    boxplot_ind(X,yl,wdth,col)
end

set(gca,'color','none','FontSize',16);
set(gca,'YTick',[],'View',[90 -90])
xlabel('R-adj')

xlim([0 1])


nbSigSess(:,2:end)./nbTotSess(:,2:end)
nbSigSess(:,2:end)

subpl=[3 4 5; 6 7 8];
for m = 1 : 2
    subplot(1,8,subpl(m,:))
    clear converge r_adj  keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    if m == 1
        col = colors(3:4,:);
    else
        col = colors(9:10,:);
    end
    xx = 0;
    for i = 1 : length(ALL)
        if  ALL(i).converge==1
            [allpb_sorted , idx] = sortrows( ALL(i).allpb(:) );
            plot(allpb_sorted,ALL(i).newf(idx),'-','Color',col(1,:));hold on
            xx = xx + 1;
            newf_avg(:,xx) = ALL(i).newf(idx);
        end
    end
    plot(allpb_sorted,mean(newf_avg'),'LineWidth',2,'Color',col(2,:))
    line([-2.5 2.5],[.5 .5],'Color','k');hold on;box on
    line([0 0],[0 1],'Color','k')

    if m == 1 
        ylabel('Probability of choosing J1')
        xlabel('log(ProbaJ1/ProbaJ2)')
    end
    set(gca,'color','none','FontSize',16);
end

% subplot(1,8,3:8)
%- plot the estimates
% nbSigSess = [];
% for m = 1 : 2
%     clear converge tStat pVal Estimate r_adj LogLik nTr keep 
%     eval(['ALL = M' num2str(m) '.ALL;']);
%     for i = 1 : length(ALL)
%         tStat(:,i) = ALL(i).mdl.Coefficients.tStat;
%         pVal(:,i) = ALL(i).mdl.Coefficients.pValue;
%         Estimate(:,i) = ALL(i).mdl.Coefficients.Estimate;
%         r_adj(:,i) = ALL(i).mdl.Rsquared.Adjusted;
%         LogLik(:,i) = ALL(i).mdl.LogLikelihood;
%         converge(:,i) = ALL(i).converge;
%         nTr(:,i) = height(ALL(i).T);
%     end
%     
%     keep = converge==1  ;
% 
%     pVal = pVal(:,keep);
%     for pp = 1 : length(pVal(:,1))
%         [h_sig, crit_p, adj_p]=fdr_bh(pVal(pp,:),0.05,'pdep','yes');
%         nbSigSess(m,pp) = sum(h_sig) ;
%         nbTotSess(m,pp) = length(h_sig) ;
%     end
% 
%     if m == 1
%         col = colors(3:4,:);
%         mm = -0.2;
%         line([0 0],[0 6],'Color','k');hold on
%     else
%         col = colors(9:10,:);
%         mm = .2;
%     end
%     for i  = 1 : length( ALL(1).mdl.CoefficientNames)-1
%         X= Estimate(i+1,keep);
%         yl=i+mm;
%         wdth = .35;
%         boxplot_ind(X,yl,wdth,col)
%     end
%     set(gca,'view',[90 -90],'color','none','FontSize',16);
%     set(gca,'YTick',1:length(ALL(1).mdl.CoefficientNames)-1,'YTickLabel',ALL(1).mdl.CoefficientNames(2:end),'YTickLabelRotation',25)
%     disp([sum(keep) length(keep)])
%     if m == 1
%         text(19.5,4,['mk M'],'FontSize',16,'Color',col(2,:))
%     else
%         text(18,4,['mk X'],'FontSize',16,'Color',col(2,:))
%     end
% end
% xlabel('Estimates')
% xlim([-20 20])
% ylim([0 length( ALL(1).mdl.CoefficientNames)])
