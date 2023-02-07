%this is a program that generates trial averages and forms clusters based
%on that. Data in TimeAnalysis/XXCBDLXXDLXX
%Created by Haleh Fotowat.

close all
clearvars
path=uigetdir;
cd(path)

%read the pooled clustering data
% load('LXLXCXCX_HA.mat')
% clear IDX
%determine what was the stimulus in between
load('Peaks_EXPFits_CBDLXXDLXX_v2.mat')
nframes=length(time);
filenames=dir('*ha_allcells.mat');
dim=NaN(1,11);
F_dff_pool_all=[];
fn=[];
stim_type_cnt=4;
dim=[];
for i=1:length(filenames)
    load(filenames(i).name,'Stim_type','Motor','F_dff_pool')
    Motor_all{i}=Motor(1:nframes);
    F_dff_pool_all=[F_dff_pool_all;F_dff_pool(:,1:nframes)];
    fn=[fn,i*ones(1,size(F_dff_pool,1))];
    if strmatch(Stim_type{21},'TM')
        dim=[dim,zeros(1,size(F_dff_pool,1))];
    else
        dim=[dim,ones(1,size(F_dff_pool,1))];

    end
end


Fit_Exp=expfit_signed;
Fit_Exp_all=expfit_signed_bf;
all_exp=[Fit_Exp(:,1);Fit_Exp(:,2);Fit_Exp(:,3);Fit_Exp(:,4)];

dfittool(all_exp)
[N,EDGES,BIN] = histcounts(all_exp,500);
NT= (EDGES(2:end)+EDGES(1:end-1))/2;
figure;bar(NT,N)
[fitresult, gof] = create2GaussFit(NT, N)

B=[fitresult.b1,fitresult.b2]
C=[fitresult.c1,fitresult.c2]
[B_s,I_s]=sort(B)
C_s=C(I_s)
%Negative threshold
th_n=B_s(2)-3*C_s(2)/(sqrt(2));
th_p=B_s(2)+3*C_s(2)/(sqrt(2));
line([th_p,th_p],[0,1800])
line([th_n,th_n],[0,1800])
% xlim([-0.06,0.06]);
% ylim([0,70])
id=1:4;
for i=1:size(F_dff_pool_all,1)
    for j=1:length(id)
        if Fit_Exp_all(i,id(j))>=th_p
            dyn_code(i,j)=2;
        elseif Fit_Exp_all(i,id(j))<=th_n
            dyn_code(i,j)=1;
        elseif ~isnan(Fit_Exp_all(i,id(j)))
             dyn_code(i,j)=0;
        else
            dyn_code(i,j)=3;
        end
    end

end

[CCStimMot,trial_ave_mat,stim_mat]=ClusterTypes(repse,stim_type_cnt,F_dff_pool_all,ISI,Motor_all,fn,fps,per2p,stim_size_degrees)

CCStimMot_th=CCStimMot(:,1:4); %do not include motor information. Failed tail tracking in all experiments of two eyes
for i=1:size(CCStimMot,1)
    for j=1:size(CCStimMot,2)-1
            if (CCStimMot(i,j))<0.2          
                CCStimMot_th(i,j)=0;
            else
                CCStimMot_th(i,j)=1;

            end
    end
    X(i)=bi2de((CCStimMot_th(i,1:4)),'left-msb');
end

[B,I]=sort(X);
Xu=unique(X);
dyn_code_m=dyn_code;
for i=1:size(dyn_code_m,1)
%     if ismember(3,dyn_code_m(i,:))%if after removing non-responsive stimuli still have no value for exp set the whole thing to 3
%         dyn_code_m(i,:)=[3,3,3,3];
%     end
    temp=num2str(dyn_code_m(i,:));
    temp=temp(find(~isspace(temp)));
    S(i)=base2dec(temp,4);
end

%% compare trials with DM to those with TM in between
%LS=11XXX
%DS=01XXX
%find the index number for LS and DS neurons under the two different
%conditions
LSD_Dind=[];
LSD_Tind=[];
DSD_Dind=[];
DSD_Tind=[];
DSI_Dind=[];
DSI_Tind=[];

for i=1:size(CCStimMot,1)
    if CCStimMot_th(i,1)==1 && CCStimMot_th(i,2)==1 && CCStimMot_th(i,3)==0  && dim(i)==1 && dyn_code(i,1)==1 && dyn_code(i,2)==1
        LSD_Dind=[LSD_Dind,i];
    elseif CCStimMot_th(i,1)==1 && CCStimMot_th(i,2)==1 && dim(i)==0 && dyn_code(i,1)==1 && dyn_code(i,2)==1
        LSD_Tind=[LSD_Tind,i];
    elseif CCStimMot_th(i,1)==0 && CCStimMot_th(i,2)==1  && dim(i)==1 && dyn_code(i,2)==1
        DSD_Dind=[DSD_Dind,i];
    elseif CCStimMot_th(i,1)==0 && CCStimMot_th(i,2)==1 && dim(i)==0 && dyn_code(i,2)==1
        DSD_Tind=[DSD_Tind,i];
    elseif CCStimMot_th(i,1)==0 && CCStimMot_th(i,2)==1 &&  dim(i)==1 && dyn_code(i,2)==2
        DSI_Dind=[DSI_Dind,i];
    elseif CCStimMot_th(i,1)==0 && CCStimMot_th(i,2)==1 && dim(i)==0 && dyn_code(i,2)==2 
        DSI_Tind=[DSI_Tind,i];
        
    end
end

figure;

hold on
title('DSI')
plot(mean(F_dff_pool_all(DSI_Dind,:)))
plot(mean(F_dff_pool_all(DSI_Tind,:)))

figure;
hold on
title('DSD')
plot(mean(F_dff_pool_all(DSD_Dind,:)))
plot(mean(F_dff_pool_all(DSD_Tind,:)))

figure;
hold on
title('LSD')
plot(mean(F_dff_pool_all(LSD_Dind,:)))
plot(mean(F_dff_pool_all(LSD_Tind,:)))

%compare the peak amplitudes
Peak_LSD=[PeakAmps(LSD_Dind,2,1);PeakAmps(LSD_Dind,2,10);PeakAmps(LSD_Dind,4,1)];
Peak_LSD_G=[ones(size(PeakAmps(LSD_Dind,2,1)));2*ones(size(PeakAmps(LSD_Dind,2,10)));3*ones(size(PeakAmps(LSD_Dind,4,1)))];
[P,A, STATS]=kruskalwallis(Peak_LSD,Peak_LSD_G)
comp_LSD_D=multcompare(STATS)
%plot this differently:
% figure
% plot(Peak_LSD_G,Peak_LSD,'*')
% hold on
% plot([1,2,3],[mean(PeakAmps(LSD_Dind,2,1)),mean(PeakAmps(LSD_Dind,2,10)),mean(PeakAmps(LSD_Dind,4,1))],'o-','LineWidth',2)
% 
PA_mean_LSD=[mean(PeakAmps(LSD_Dind,2,1)),mean(PeakAmps(LSD_Dind,2,10)),mean(PeakAmps(LSD_Dind,4,1))];
PA_std_LSD=[std(PeakAmps(LSD_Dind,2,1))/sqrt(length(PeakAmps(LSD_Dind,2,1))),std(PeakAmps(LSD_Dind,2,10))/sqrt(length(PeakAmps(LSD_Dind,2,10))),std(PeakAmps(LSD_Dind,4,1))/sqrt(length(PeakAmps(LSD_Dind,4,1)))];
figure
bar(PA_mean_LSD)
hold on
errorbar(PA_mean_LSD,PA_std_LSD,'.')

figure
subplot(3,1,1)
histogram(PeakAmps(LSD_Dind,2,1),'Normalization','pdf')
subplot(3,1,2)
histogram(PeakAmps(LSD_Dind,2,10),'Normalization','pdf')
subplot(3,1,3)
histogram(PeakAmps(LSD_Dind,4,1),'Normalization','pdf')


index=find(PeakAmps(LSD_Dind,4,1)<0.01);
figure;plot(mean(F_dff_pool_all(LSD_Dind(index),:)))

index=find(PeakAmps(LSD_Tind,4,1)<0.01);
figure;plot(mean(F_dff_pool_all(LSD_Tind(index),:)))


Peak_LSD=[PeakAmps(LSD_Tind,2,1);PeakAmps(LSD_Tind,2,10);PeakAmps(LSD_Tind,4,1)];
Peak_LSD_G=[ones(size(PeakAmps(LSD_Tind,2,1)));2*ones(size(PeakAmps(LSD_Tind,2,10)));3*ones(size(PeakAmps(LSD_Tind,4,1)))];
[P,A, STATS]=kruskalwallis(Peak_LSD,Peak_LSD_G)

% figure
% plot(Peak_LSD_G,Peak_LSD,'*')
% hold on
% plot([1,2,3],[mean(PeakAmps(LSD_Tind,2,1)),mean(PeakAmps(LSD_Tind,2,10)),mean(PeakAmps(LSD_Tind,4,1))],'o-','LineWidth',2)

PA_mean_LST=[mean(PeakAmps(LSD_Tind,2,1));mean(PeakAmps(LSD_Tind,2,10));mean(PeakAmps(LSD_Tind,4,1))]
PA_std_LST=[std(PeakAmps(LSD_Tind,2,1))/sqrt(length(PeakAmps(LSD_Tind,2,1)));std(PeakAmps(LSD_Tind,2,10))/sqrt(length(PeakAmps(LSD_Tind,2,1)));std(PeakAmps(LSD_Tind,4,1))/sqrt(length(PeakAmps(LSD_Tind,2,1)))]
figure
bar(PA_mean_LST)
hold on
errorbar(PA_mean_LST,PA_std_LST,'.')

%% bar plot LST and LSD together in one 
figure
bar([PA_mean_LST',PA_mean_LSD])
hold on
errorbar([PA_mean_LST',PA_mean_LSD],[PA_std_LST',PA_std_LSD],'.')
%compare trial 10 and 11 for dimming
Peak_LSD=[PeakAmps(LSD_Tind,2,10);PeakAmps(LSD_Tind,4,1)];
Peak_LSD_G=[ones(size(PeakAmps(LSD_Tind,2,1)));2*ones(size(PeakAmps(LSD_Tind,4,1)))];
[P,A, STATS]=kruskalwallis(Peak_LSD,Peak_LSD_G)
Peak_LST=[PeakAmps(LSD_Dind,2,10);PeakAmps(LSD_Dind,4,1)];
Peak_LST_G=[ones(size(PeakAmps(LSD_Dind,2,1)));2*ones(size(PeakAmps(LSD_Dind,4,1)))];
[P,A, STATS]=kruskalwallis(Peak_LST,Peak_LST_G)

Peak_LSTD=[PeakAmps(LSD_Tind,4,1);PeakAmps(LSD_Dind,4,1)];
Peak_LSTD_G=[ones(size(PeakAmps(LSD_Tind,4,1)));2*ones(size(PeakAmps(LSD_Dind,4,1)))];
[P,A, STATS]=kruskalwallis(Peak_LSTD,Peak_LSTD_G)

Peak_LSTD=[PeakAmps(LSD_Tind,1,1);PeakAmps(LSD_Dind,1,1)];
Peak_LSTD_G=[ones(size(PeakAmps(LSD_Tind,1,1)));2*ones(size(PeakAmps(LSD_Dind,1,1)))];
[P,A, STATS]=kruskalwallis(Peak_LSTD,Peak_LSTD_G)

%%

figure
subplot(3,1,1)
histogram(PeakAmps(LSD_Tind,2,1),'Normalization','pdf')
subplot(3,1,2)
histogram(PeakAmps(LSD_Tind,2,10),'Normalization','pdf')
subplot(3,1,3)
histogram(PeakAmps(LSD_Tind,4,1),'Normalization','pdf')

comp_LSD_T=multcompare(STATS)


%% plot DSI during time in between
figure;
stdshade(F_dff_pool_all(DSI_Tind,fps*10+1:end),0.2,'k',time(1:fps*30))
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','w','r'};
hold on
for k=1:3
    %plot repse areas at first color
    for kk=1:10
        if kkk==1
            extra_time=0;
        else
            extra_time=5*per2p;
        end
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight

%% plot LSD time course for the figure
figure;
subplot(2,1,1)
stdshade(F_dff_pool_all(LSD_Tind,fps*10+1:end),0.2,'r',time(1:fps*30))
hold on
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','y','r'};

for k=1:3
    %plot repse areas at first color
    for kk=1:10
        if kkk==1
            extra_time=0;
        else
            extra_time=5*per2p;
        end
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight
subplot(2,1,2)
hold on;stdshade(F_dff_pool_all(LSD_Dind,fps*10+1:end),0.2,'k',time(1:fps*30))
hold on
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','k','r'};

for k=1:3
    %plot repse areas at first color
    for kk=1:10
        if kkk==1
            extra_time=0;
        else
            extra_time=5*per2p;
        end
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight



%% not in the analysis



% figure;plot(mean(F_dff_pool_all(CCStimMot_th(:,1)==1 & CCStimMot_th(:,2)==1  & CCStimMot_th(:,3)==0 & dim'==1,:)))
% hold on
% plot(mean(F_dff_pool_all(CCStimMot_th(:,1)==1 & CCStimMot_th(:,2)==1  & CCStimMot_th(:,3)==1 & dim'==1,:)))

%use kmeans to divide the clusters in half to see if I can detect
%everything responsive neurons.
% indexD=find(CCStimMot_th(:,1)==1 & CCStimMot_th(:,2)==1 & dim'==1 & dyn_code(:,1)==1 & dyn_code(:,2)==1);%pick all the looming responsive decreasing neurons
indexD=find(CCStimMot_th(:,1)==1 & CCStimMot_th(:,2)==1 & dim'==1);%pick all the looming responsive decreasing neurons

PTimes=nanmean(PeakIndrelstim(indexD,2,:),3);%run k means to separate two clusters based on peak timing. This is based on my prior knowledge that the looming-specific neurons peak earlier than non-specific ones
IDXD=kmeans(PTimes,2);
%compare peak times of the two categories 
[P,A,STAT]=kruskalwallis(PTimes,IDXD)
[a,minin]=min(STAT.meanranks);
[a,maxin]=max(STAT.meanranks);
%assign the pure looming-sensitive decreasing neurons to the cluster with
%the earlier peak. and the non-specific looming sensitive decreasing
%neurons to the other cluster.
PLSD_D_ind=indexD(IDXD==minin);%pure LSD dim index
NLSD_D_ind=indexD(IDXD==maxin);%pure LSD dim index


%Now use the same method to isolate the pure looming sensitive neurons form
%none-specific ones even when the dimming stimulus is not present (anyway
%the selection is based on the peak time for DL stimulus and does not
%depend on knowing anything about the dimming response.
%indexT=find(CCStimMot_th(:,1)==1 & CCStimMot_th(:,2)==1 & dim'==0 & (dyn_code(:,1)==1 | dyn_code(:,2)==1));
indexT=find(CCStimMot_th(:,1)==1 & CCStimMot_th(:,2)==1 & dim'==0);

PTimes=nanmean(PeakIndrelstim(indexT,2,:),3);
IDXT=kmeans(PTimes,2);
%compare peak times of the two categories
[P,A,STAT]=kruskalwallis(PTimes,IDXT)
[a,minin]=min(STAT.meanranks);
[a,maxin]=max(STAT.meanranks);

PLSD_T_ind=indexT(IDXT==minin);%pure LSD dim index
NLSD_T_ind=indexT(IDXT==maxin);%pure LSD dim index

figure;stdshade(F_dff_pool_all(PLSD_T_ind,:),0.2,'k',time)
hold on;stdshade(F_dff_pool_all(PLSD_D_ind,:),0.2,'r',time)
hold on
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'m','r','c','r'};

for k=1:4
    %plot repse areas at first color
    for kk=1:10
        if kkk==1
            extra_time=0;
        else
            extra_time=5*per2p;
        end
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight
length(PLSD_T_ind)
length(PLSD_D_ind)
% figure;plot(mean(F_dff_pool_all(PLSD_D_ind,:)))
% hold on;plot(mean(F_dff_pool_all(NLSD_D_ind,:)))

%% form ratio of the response in the looming sensitive neurons
Ratio_PLSD_D=(PeakAmps(PLSD_D_ind,4,1)-PeakAmps(PLSD_D_ind,2,1))./(PeakAmps(PLSD_D_ind,4,1)+PeakAmps(PLSD_D_ind,2,1));
Ratio_PLSD_T=(PeakAmps(PLSD_T_ind,4,1)-PeakAmps(PLSD_T_ind,2,1))./(PeakAmps(PLSD_T_ind,4,1)+PeakAmps(PLSD_T_ind,2,1));

Ratios_PLSD=[Ratio_PLSD_D;Ratio_PLSD_T];
Groups_PLSD=[ones(size(Ratio_PLSD_D));2*ones(size(Ratio_PLSD_T))];
p=kruskalwallis(Ratios_PLSD,Groups_PLSD)
title(strcat('p=' ,num2str(p),', Nd(1)=',num2str(length(PLSD_D_ind)),', Nt(2)=',num2str(length(PLSD_T_ind))))

%% form ratio of the response in the dimming sensitive decreasing neurons
%see if these can be divided into two groups those that respond both to DL
%and DM and those that respond only to DL . It actually did not show a big
%difference in tuning. They all responded both to dimming and dark looming
%so merging the two
% PTimes=nanmean(PeakIndrelstim(DSD_Dind,2,:),3);%run k means to separate two clusters based on peak timing. This is based on my prior knowledge that the looming-specific neurons peak earlier than non-specific ones
% IDXD=kmeans(PTimes,2);
% %compare peak times of the two categories 
% [P,A,STAT]=kruskalwallis(PTimes,IDXD)
% [a,minin]=min(STAT.meanranks);
% [a,maxin]=max(STAT.meanranks);
% %assign the pure looming-sensitive decreasing neurons to the cluster with
% %the earlier peak. and the non-specific looming sensitive decreasing
% %neurons to the other cluster.
% PDSD_D_ind=DSD_Dind(IDXD==minin);%pure LSD dim index
% NDSD_D_ind=DSD_Dind(IDXD==maxin);%pure LSD dim index
% figure;stdshade(F_dff_pool_all(PDSD_D_ind,:),0.2,'k',time)
% hold on;stdshade(F_dff_pool_all(NDSD_D_ind,:),0.2,'r',time)
% 
% 


Ratio_DSD_D=(PeakAmps(DSD_Dind,4,1)-PeakAmps(DSD_Dind,2,1))./(PeakAmps(DSD_Dind,4,1)+PeakAmps(DSD_Dind,2,1));
Ratio_DSD_T=(PeakAmps(DSD_Tind,4,1)-PeakAmps(DSD_Tind,2,1))./(PeakAmps(DSD_Tind,4,1)+PeakAmps(DSD_Tind,2,1));

Ratios_DSD=[Ratio_DSD_D;Ratio_DSD_T];
Groups_DSD=[ones(size(Ratio_DSD_D));2*ones(size(Ratio_DSD_T))];
p=kruskalwallis(Ratios_DSD,Groups_DSD)
title(strcat('p=' ,num2str(p),', Nd(1)=',num2str(length(DSD_Dind)),', Nt(2)=',num2str(length(DSD_Tind))))

%% form ratio of the response in the dimming sensitive increasing neurons

Ratio_DSI_D=(PeakAmps(DSI_Dind,4,1)-PeakAmps(DSI_Dind,2,1))./(PeakAmps(DSI_Dind,4,1)+PeakAmps(DSI_Dind,2,1));
Ratio_DSI_T=(PeakAmps(DSI_Tind,4,1)-PeakAmps(DSI_Tind,2,1))./(PeakAmps(DSI_Tind,4,1)+PeakAmps(DSI_Tind,2,1));

Ratios_DSI=[Ratio_DSI_D;Ratio_DSI_T];
Groups_DSI=[ones(size(Ratio_DSI_D));2*ones(size(Ratio_DSI_T))];
p=kruskalwallis(Ratios_DSI,Groups_DSI)
title(strcat('p=' ,num2str(p),', Nd(1)=',num2str(length(DSI_Dind)),', Nt(2)=',num2str(length(DSI_Tind))))


% %% Plot the trial average clusters
% figure;
% subplot(1,6,1)
% %DM
% ax1=imagesc(trial_ave_mat(I,1:(fps+1)))
% hold on
% line([f0,f0],[1,length(I)],'color','w')
% line([(fps+1)-f0,(fps+1)-f0],[1,length(I)],'color','w')
% %colorbar
% caxis([-2,6])
% 
% subplot(1,6,2)
% %DL
% ax2=imagesc(trial_ave_mat(I,(fps+1):(fps+1)*2))
% line([f0,f0],[1,length(I)],'color','w')
% line([(fps+1)-f0,(fps+1)-f0],[1,length(I)],'color','w')
% %colorbar
% caxis([-2,6])
% 
% 
% subplot(1,6,3)
% %CB
% ax3=imagesc(trial_ave_mat(I,(fps+1)*2:(fps+1)*3)) %start 10 frames before stim onset
% hold on
% line([f0,f0],[1,length(I)],'color','w')
% line([(fps+1)-f0,(fps+1)-f0],[1,length(I)],'color','w')
% % colorbar
% caxis([-2,6])
% 
% 
% subplot(1,6,4)
% %CB
% ax4=imagesc(trial_ave_mat(I,(fps+1)*3:end)) %start 10 frames before stim onset
% hold on
% line([f0,f0],[1,length(I)],'color','w')
% line([(fps+1)-f0,(fps+1)-f0],[1,length(I)],'color','w')
% %colorbar
% caxis([-2,6])
% 
% subplot(1,6,5)
% %CB
% ax5=imagesc(CCStimMot_th(I,:)) %start 10 frames before stim onset
% hold on
% 
% ax6=subplot(1,6,6)
% %CB
% imagesc(X(I)');colormap(ax6,1-jet) %start 10 frames before stim onset
% hold on
% 
% %plot the time course
% figure
% index=1:100:length(X);
% colormap1=1-jet;
% 
% for i=1:length(index)-1
%     %figure
%     subplot(length(index)-1,1,i)
%     hold on
%     plot(1:(fps+1),median(trial_ave_mat(I(index(i):index(i+1)-1),1:(fps+1)),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:))
%     line([f0,f0],[-0.5,0.5],'color','k')
%     line([(fps+1)-f0,(fps+1)-f0],[-0.5,0.5],'color','k')
%     plot((fps+1)+1:(fps+1)*2,median(trial_ave_mat(I(index(i):index(i+1)-1),(fps+1)+1:(fps+1)*2),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:))
%     line([(fps+1)+f0,(fps+1)+f0],[-0.5,0.5],'color','k')
%     line([(fps+1)*2-f0,(fps+1)*2-f0],[-0.5,0.5],'color','k')
% 
%     plot((fps+1)*2+1:(fps+1)*3,median(trial_ave_mat(I(index(i):index(i+1)-1),(fps+1)*2+1:(fps+1)*3),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:))
%     line([(fps+1)*2+f0,(fps+1)*2+f0],[-0.5,0.5],'color','k')
%     line([(fps+1)*3-f0,(fps+1)*3-f0],[-0.5,0.5],'color','k')
%     
%     plot((fps+1)*3+1:(fps+1)*4,median(trial_ave_mat(I(index(i):index(i+1)-1),(fps+1)*3+1:end),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:))
%     line([(fps+1)*3+f0,(fps+1)*3+f0],[-0.5,0.5],'color','k')
%     line([(fps+1)*4-f0,(fps+1)*4-f0],[-0.5,0.5],'color','k')
%     %axis tight
%     xlim([1,(fps+1)*4])
%     ylim([-1,3])
%     axis off
%     
%     
% end
% 
% %define various response categories
% Xu=unique(X);
% for i=1:length(Xu)
%     XuL(i)=sum(X==Xu(i));
%     fnX(i)=length(unique(fn(X==Xu(i))));
% end
% %take only units that were present in at least 3 fish
% 
% %take only categories that contained more than 100 neurons
% Index=find(XuL>50);
% for i=2:length(Index)
%     figure;
%     
%     title(strcat('Clust=',num2str(Xu(Index(i))),', n=',num2str(sum(X==Xu(Index(i))))))
%     hold on
%     stdshade(F_dff_pool_all(X==Xu(Index(i)),:),0.2,'k',time)
%     hold on
%     LIM=ylim(gca);
%     DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
%     kkk=1;
%     color_bar={'r','m','b','c'};
%     
%     for k=1:4
%         %plot repse areas at first color
%         for kk=1:10
%             if kkk==1
%                 extra_time=0;
%             else
%                 extra_time=5*per2p;
%             end
%             area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
%             area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
%             
%             kkk=kkk+1;
%         end
%         
%     end
%     axis tight
% end
% colormap1=1-jet;
% for i=2:length(Index)
%     figure;
%     title(strcat('Clust=',num2str(Xu(Index(i))),', n=',num2str(sum(X==Xu(Index(i))))))
%     color=colormap1(Xu(Index(i))*17+1,:);
%     hold on
%     plot(time(f0:f0+fps),median(trial_ave_mat(X==Xu(Index(i)),1:(fps+1)),1),'Color',colormap1(17*Xu(Index(i))+1,:))
%     hold on
%     plot(time(f0+(fps+1):f0+(fps+1)*2),median(trial_ave_mat(X==Xu(Index(i)),fps+1:(fps+1)*2),1),'Color',colormap1(17*Xu(Index(i))+1,:))
%     plot(time(f0+(fps+1)*2+1:f0+(fps+1)*3),median(trial_ave_mat(X==Xu(Index(i)),(fps+1)*2+1:(fps+1)*3),1),'Color',colormap1(17*Xu(Index(i))+1,:))
%     plot(time(f0+(fps+1)*3+1:f0+(fps+1)*4),median(trial_ave_mat(X==Xu(Index(i)),(fps+1)*3+1:end),1),'Color',colormap1(17*Xu(Index(i))+1,:))
%     axis tight
%     LIM=ylim(gca);
%     DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
%     kkk=1;
%     color_bar={'k','r','c','r','b'};
%     
%     for k=1:4
%         %plot repse areas at first color
%         for kk=1:1
%           
%             extra_time=5*per2p;
%             area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
%             area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
%             
%             kkk=kkk+1;
%         end
%         
%     end
% end
% PeakInd_relon=NaN(size(PeakInd));
% %compare response timing for the trial ave across different 
% filenames=dir('2*AllBrains_procdata_ha_allcells.mat');
% 
% for i=1:length(filenames)
%     load(filenames(i).name,'stim_size_degrees');
%     theta_i=min(stim_size_degrees)/2;
%     TTC(i)=-1*480/tand(theta_i);
% end
% for k=1:size(trial_ave_mat,1)
%     for i=1:stim_type_cnt
%         lb=1+(fps+1)*(i-1);
%         rb=(fps+1)*i;
%         [m,in]=findpeaks(trial_ave_mat(k,lb:rb));
%         if ~isempty(m)
%             [a,b]=max(m);
%             trial_ave_pt(k,i)=in(b)-f0;
%             trial_ave_ttc(k,i)=trial_ave_pt(k,i)*per2p+TTC(fn(k))/1000;
%         else
%             trial_ave_pt(k,i)=NaN;
%         end
%     end
% end
% 
% Peak_times=[trial_ave_ttc(X==12,2);trial_ave_ttc(X==5,2);trial_ave_ttc(X==10,2)];
% groups=[ones(size(trial_ave_ttc(X==12,2)));2*ones(size(trial_ave_ttc(X==5,2)));3*ones(size(trial_ave_ttc(X==10,2)))];
% [p,a,STATS]=kruskalwallis(Peak_times,groups);
% multcompare(STATS)


IDX=X;
save('Peaks_EXPFits_ClustsNew_CBDLXXDLXX_v2.mat')

% hold on
% LIM=ylim(gca);
% DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
% kkk=1;
% color_bar={'k','r','c','r','b'};
% 
% for k=1:4
%     %plot repse areas at first color
%     for kk=1:10
%         
%         extra_time=5*per2p;
%         area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
%         area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
%         
%         kkk=kkk+1;
%     end
%     
% end
