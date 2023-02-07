%this is a program that generates trial averages and forms clusters based
%on that. Created by Haleh Fotowat.
close all
clearvars
path=uigetdir;
cd(path)

%read the pooled clustering data
load('DMDLBRDLCB_082628_091112_HA.mat')
clear IDX
load('Peaks_EXPFits_DMDLBRDLCB_v2.mat')
ISI=40;
per2p=1.06295;
stim_type_cnt=5;
repse=5;
fps=69;
f0=round(ISI/(2*per2p));%half of tthe ISI
nframes=size(F_dff_pool_all,2)
framerate=1/per2p;

time=1/framerate:1/framerate:nframes/framerate;


[CCStimMot,trial_ave_mat,stim_mat]=ClusterTypes(repse,stim_type_cnt,F_dff_pool_all,ISI,Motor_all,fn,fps,per2p,stim_size_degrees)
CCStimMot_th=CCStimMot;
for i=1:size(CCStimMot,1)
    for j=1:size(CCStimMot,2)
        if j~=6
            if abs(CCStimMot(i,j))<0.2          
                CCStimMot_th(i,j)=0;
            else
                CCStimMot_th(i,j)=1;

            end
        else
            if (CCStimMot(i,j))<0.5         
                CCStimMot_th(i,j)=0;
            else
                CCStimMot_th(i,j)=1;

            end
        end
    end
    X(i)=bi2de((CCStimMot_th(i,[1:3,5])),'left-msb');
    
end
mask=CCStimMot_th(:,[1:2,4:5]);


[B,I]=sort(X);
Xu=unique(X);
IDX=X;




%% check and see if the timing of response changes over trials for different looming-responsive clusters
for i=1:size(F_dff_pool_all,1)
    PeakTime_DL(i,:)=squeeze(PeakIndrelstim(i,2,:));
end
PeakTime_LS=[PeakTime_DL(X==5,2);PeakTime_DL(X==5,3);PeakTime_DL(X==5,4)];
groups=[2*ones(size(PeakTime_DL(X==5,2)));3*ones(size(PeakTime_DL(X==5,3)));4*ones(size(PeakTime_DL(X==5,4)))];
[a,b,STAS]=kruskalwallis(PeakTime_LS,groups)
%% Plot the trial average clusters
figure;
subplot(1,6,1)
%DM
ax1=imagesc(trial_ave_mat(I,1:70))
hold on
line([f0,f0],[1,length(I)],'color','w')
line([70-f0,70-f0],[1,length(I)],'color','w')
%colorbar
caxis([-2,6])

subplot(1,6,2)
%DL
ax2=imagesc(trial_ave_mat(I,70:70*2))
line([f0,f0],[1,length(I)],'color','w')
line([70-f0,70-f0],[1,length(I)],'color','w')
%colorbar
caxis([-2,6])


subplot(1,6,3)
%CB
ax3=imagesc(trial_ave_mat(I,70*2:70*3)) %start 10 frames before stim onset
hold on
line([f0,f0],[1,length(I)],'color','w')
line([70-f0,70-f0],[1,length(I)],'color','w')
% colorbar
caxis([-2,6])


subplot(1,6,4)
%CB
ax4=imagesc(trial_ave_mat(I,70*4:end)) %start 10 frames before stim onset
hold on
line([f0,f0],[1,length(I)],'color','w')
line([70-f0,70-f0],[1,length(I)],'color','w')
%colorbar
caxis([-2,6])

% subplot(1,7,5)
% caxis([-2,6])
% colorbar

subplot(1,6,5)
%CB
ax5=imagesc(CCStimMot_th(I,[1:3,5])) %start 10 frames before stim onset
hold on

ax6=subplot(1,6,6)
%CB
imagesc(X(I)');colormap(ax6,1-jet) %start 10 frames before stim onset
hold on


figure;
caxis([-2,6])
colorbar
figure;
[B,Im]=sort(CCStimMot_th(:,6));
All_sens=(CCStimMot_th(:,1)==1 | CCStimMot_th(:,2)==1 | CCStimMot_th(:,3)==1 | CCStimMot_th(:,5)==1);
Sensory_Motor=[All_sens(Im),CCStimMot_th(Im, 6)];
imagesc(Sensory_Motor) %start 10 frames before stim onset


%% plot the time course of binned clusters
figure
index=1:100:length(X);
colormap1=1-jet;
for i=1:length(index)-1
    %figure
    subplot(length(index)-1,1,i)
    hold on
    plot(1:70,median(trial_ave_mat(I(index(i):index(i+1)-1),1:70),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:),'LineWidth',1.5)
    line([f0,f0],[-0.5,0.5],'color','k')
    line([70-f0,70-f0],[-0.5,0.5],'color','k')
    plot(71:140,median(trial_ave_mat(I(index(i):index(i+1)-1),70+1:70*2),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:),'LineWidth',1.5)
    line([70+f0,70+f0],[-0.5,0.5],'color','k')
    line([140-f0,140-f0],[-0.5,0.5],'color','k')

    plot(141:210,median(trial_ave_mat(I(index(i):index(i+1)-1),70*2+1:70*3),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:),'LineWidth',1.5)
    line([140+f0,140+f0],[-0.5,0.5],'color','k')
    line([210-f0,210-f0],[-0.5,0.5],'color','k')
    
    plot(211:280,median(trial_ave_mat(I(index(i):index(i+1)-1),70*4+1:end),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:),'LineWidth',1.5)
    line([211+f0,211+f0],[-0.5,0.5],'color','k')
    line([280-f0,280-f0],[-0.5,0.5],'color','k')
    %axis tight
    xlim([1,280])
    ylim([-2,3])
    axis off
    
    
end
%define various response categories
for i=1:length(Xu)
    XuL(i)=sum(X==Xu(i));
    fnX(i)=length(unique(fn(X==Xu(i))));
end

Index=find(XuL>1);
colormap1=1-jet;
for i=1:length(Index)
    figure;
    title(strcat('Clust=',num2str(Xu(Index(i))),', n=',num2str(sum(X==Xu(Index(i))))))
    color=colormap1(Xu(Index(i))*17+1,:);
    hold on
    %stdshade((F_dff_pool_all(X==Xu(Index(i)),:)),0.2,'k',time)
    stdshade((F_dff_pool_all(X==Xu(Index(i)),[1:15*fps,20*fps+1:end])),0.2,'k',time(1:20*fps))

    hold on
    %ylim([-0.5,5])
    LIM=ylim(gca);
    DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
    kkk=1;
    %color_bar={'k','r','c','r','b'};
    color_bar={'k','r','c','b'};
   
    for k=1:4%5
        %plot repse areas at first color
        for kk=1:5
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
    
end

colormap1=1-jet;
for i=1:length(Index)
    figure;
    title(strcat('Clust=',num2str(Xu(Index(i))),', n=',num2str(sum(X==Xu(Index(i))))))
    color=colormap1(Xu(Index(i))*17+1,:);
    hold on
    plot(time(f0:f0+69),median(trial_ave_mat(X==Xu(Index(i)),1:70),1),'Color',colormap1(17*Xu(Index(i))+1,:))
    hold on
    plot(time(f0+70:f0+139),median(trial_ave_mat(X==Xu(Index(i)),71:140),1),'Color',colormap1(17*Xu(Index(i))+1,:))
    plot(time(f0+140:f0+209),median(trial_ave_mat(X==Xu(Index(i)),141:210),1),'Color',colormap1(17*Xu(Index(i))+1,:))
    plot(time(f0+210:f0+279),median(trial_ave_mat(X==Xu(Index(i)),70*4+1:end),1),'Color',colormap1(17*Xu(Index(i))+1,:))
    axis tight
    hold on
    LIM=ylim(gca);
    DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
    kkk=1;
    color_bar={'k','r','c','r','b'};
    
    for k=[1:3,5]
        %plot repse areas at first color
        for kk=1:1
          
            extra_time=5*per2p;
            area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
            area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
            
            kkk=kkk+1;
        end
        
    end
end
%% plot all units that  respond to either dimming, or looming (CB or DL)
Index=Xu(2:end);
figure; hold on
colormap1=1-jet;
for i=1:length(Index)
    figure;
    hold on
    %title(strcat('Clust=',num2str((Index(i))),', n=',num2str(sum(X==(Index(i))))))
    color=colormap1((Index(i))*17+1,:);
    hold on
    plot(time(f0:f0+69),median(trial_ave_mat(X==(Index(i)),1:70),1),'Color',colormap1(17*(Index(i))+1,:),'LineWidth',1)
    hold on
    plot(time(f0+70:f0+139),median(trial_ave_mat(X==(Index(i)),71:140),1),'Color',colormap1(17*(Index(i))+1,:),'LineWidth',1)
    plot(time(f0+140:f0+209),median(trial_ave_mat(X==(Index(i)),141:210),1),'Color',colormap1(17*(Index(i))+1,:),'LineWidth',1)
    plot(time(f0+210:f0+279),median(trial_ave_mat(X==(Index(i)),70*4+1:end),1),'Color',colormap1(17*(Index(i))+1,:),'LineWidth',1)
    axis tight
    LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=[1:3,5]
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end

title(strcat('Cluster',num2str(Index(i))))
end




%calculate time to collision at stimulus onset
filenames=dir('2*AllBrains_procdata_ha_allcells.mat');
for i=1:length(filenames)
    load(filenames(i).name,'stim_size_degrees');
    theta_i=min(stim_size_degrees)/2;
    TTC(i)=-1*480/tand(theta_i);
end
for k=1:size(trial_ave_mat,1)
    for i=1:stim_type_cnt
        lb=1+(fps+1)*(i-1);
        rb=(fps+1)*i;
        [m,in]=findpeaks(trial_ave_mat(k,lb:rb));
        if ~isempty(m)
            [a,b]=max(m);
            trial_ave_pt(k,i)=in(b)-f0;
            trial_ave_ttc(k,i)=trial_ave_pt(k,i)*per2p+TTC(fn(k))/1000;
        else
            trial_ave_pt(k,i)=NaN;
        end
    end
end

Peak_times=[trial_ave_ttc(X==12,2);trial_ave_ttc(X==5,2);trial_ave_ttc(X==4,2);trial_ave_ttc(X==3,2);trial_ave_ttc(X==1,2)];
groups=[ones(size(trial_ave_ttc(X==12,2)));2*ones(size(trial_ave_ttc(X==5,2)));3*2*ones(size(trial_ave_ttc(X==4,2)));4*2*ones(size(trial_ave_ttc(X==3,2)));5*2*ones(size(trial_ave_ttc(X==1,2)))];
[p,a,STATS]=kruskalwallis(Peak_times,groups);
figure
multcompare(STATS)
%% compare all peak times
Peak_times=[];
groups=[];
index=[5,12];
for i=1:length(index)
    Peak_times=[Peak_times;trial_ave_ttc(X==index(i),2)];
    groups=[groups;i*ones(size(trial_ave_ttc(X==index(i),2)))];
end
[p,a,STATS]=kruskalwallis(Peak_times,groups);
figure;
multcompare(STATS)

figure;
bar([5,12],[mean(trial_ave_ttc(X==5,2)),mean(trial_ave_ttc(X==12,2))])
hold on
errorbar([5,12],[mean(trial_ave_ttc(X==5,2)),mean(trial_ave_ttc(X==12,2))],[std(trial_ave_ttc(X==5,2)),std(trial_ave_ttc(X==12,2))],'.')

figure;
bar([5,12],[mean(trial_ave_pt(X==5,2)),mean(trial_ave_pt(X==12,2))])
hold on
errorbar([5,12],[mean(trial_ave_pt(X==5,2)),mean(trial_ave_pt(X==12,2))],[std(trial_ave_pt(X==5,2)),std(trial_ave_pt(X==12,2))],'.')
%% plot example average trace and bar plot for time comparison
figure
subplot(2,1,1)
plot(time(f0:f0+fps),stim_size_degrees(f0:f0+fps),'k','LineWidth',1)
hold on
axis tight
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r'};

for k=1
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        %area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
hold on
subplot(2,1,2)
stdshade(trial_ave_mat(X==5,71:140),0.2,'r',time(f0:f0+69))
hold on
stdshade(trial_ave_mat(X==12,71:140),0.2,'b',time(f0:f0+69))

axis tight
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r'};

for k=1
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end



%show percentage of inhibitor neurons.
RegData=csvread('DMDLBRDLCB_celltypes_registered_v2.csv');
EI_all=RegData(:,7);

figure;
bar([5,12],[100*sum(IDX==5& EI_all'==-1)/sum(IDX==5),100*sum(IDX==12& EI_all'==-1)/sum(IDX==12)])
hold on

%% calculate dynamic clusters
%set the exponential code for the non-relevant stimulus to zero
% dyn_code_m=mask.*dyn_code;
Fit_Exp=expfit_signed;
Fit_Exp_all=expfit_signed_bf;%includes bad fits
exp_CB=expfit_signed(CCStimMot_th(:,5)==1,5);
exp_DL=expfit_signed(CCStimMot_th(:,2)==1,2);
exp_DM=expfit_signed(CCStimMot_th(:,1)==1,1);

all_exp=[exp_CB;exp_DL;exp_DM];
dfittool(all_exp)
figure;histogram(all_exp,100,'Normalization','count');


[N,EDGES] = histcounts(all_exp,10000);

%remove the near zero bars
NT= (EDGES(2:end)+EDGES(1:end-1))/2;
rb=find(NT>0,1);
lb=rb-1;
figure;bar(NT,N)
% cftool(NT,N)
figure
[fitresult, gof] = create2GaussFit(NT, N)
%c parameter is sqrt(2)*sigma, so divide the c parameter by sqrt(2) to get
%the standard deviation and go one std away from the mean for threshold.
B=[fitresult.b1,fitresult.b2]
C=[fitresult.c1,fitresult.c2]
[B_s,I_s]=sort(B)
C_s=C(I_s)

figure;plot(NT,fitresult(NT),'r')
% hold on
% % bar(NT,N)
% breakplot(NT,fitresult(NT),2,120,'Line')
% hold on

axis tight
%Negative threshold
FWHM=2*sqrt(2*log(2))*(C_s(2)/sqrt(2))
th_n=B_s(2)-FWHM/2;
th_p=B_s(2)+FWHM/2;
line([th_p,th_p],[0,max(N)])
line([th_n,th_n],[0,max(N)])
%% Figure 6a
% figure;plot(NT(NT<th_n),N(NT<th_n),'Color',[0.4660 0.6740 0.1880],'Marker','.')
% hold on
% plot(NT(NT>=th_p),N(NT>=th_p),'Color',[0.6350 0.0780 0.1840],'Marker','.')
% plot(NT(NT>=th_n&NT<th_p),N(NT>=th_n&NT<th_p),'k.')
figure
 [fitresult, gof] = create2GaussFit(NT, N)
%xlim([LeftB,RightB])
line([th_p,th_p],[0,max(N)])
line([th_n,th_n],[0,max(N)])
hold on
plot(NT(NT<th_n),N(NT<th_n),'Color',[0.4660 0.6740 0.1880],'Marker','.','LineStyle','none')
hold on
plot(NT(NT>=th_p),N(NT>=th_p),'Color',[0.6350 0.0780 0.1840],'Marker','.','LineStyle','none')
plot(NT(NT>=th_n&NT<th_p),N(NT>=th_n&NT<th_p),'k.')
axis tight
legend off
box off
%zoomed version
figure
[fitresult, gof] = create2GaussFit(NT, N)
%xlim([LeftB,RightB])
line([th_p,th_p],[0,max(N)])
line([th_n,th_n],[0,max(N)])
hold on
plot(NT(NT<th_n),N(NT<th_n),'Color',[0.4660 0.6740 0.1880],'Marker','.','LineStyle','none')
hold on
plot(NT(NT>=th_p),N(NT>=th_p),'Color',[0.6350 0.0780 0.1840],'Marker','.','LineStyle','none')
plot(NT(NT>=th_n&NT<th_p),N(NT>=th_n&NT<th_p),'k.','MarkerSize',1)
axis tight
%xlim([-0.5e-3,0.5e-3])
xlim([-0.001,0.001])

legend off
box off

%percentage in each group
Dec=sum(N(NT<th_n))/sum(N)
Pot=sum(N(NT>=th_p))/sum(N)
Stab=sum(N(NT>=th_n&NT<th_p))/sum(N)
%%
% ylim([0,50])
% th_n0=B_s(2)-FWHM/2;
% th_p0=B_s(2)+FWHM/2;
%dfittool(all_exp)
figure;histogram(all_exp,10000,'Normalization','count')
hold on
line([th_p,th_p],[0,max(N)])
line([th_n,th_n],[0,max(N)])

[N1,EDGES1]=histcounts(all_exp,10000);
[a,b]=max(N1);
LeftB=EDGES1(b);
RightB=EDGES1(b+1);
y=[sum(all_exp>=LeftB & all_exp<th_n),sum(all_exp>=th_n & all_exp<th_p),sum(all_exp>=th_p & all_exp<RightB)]
bar(mean([LeftB,RightB]),y,RightB-LeftB, 'stacked')

xlim([-0.1,0.1])

figure
[fitresult, gof] = create2GaussFit(NT, N)
xlim([LeftB,RightB])
line([th_p,th_p],[0,max(N)])
line([th_n,th_n],[0,max(N)])

figure
[fitresult, gof] = create2GaussFit(NT, N)
xlim([LeftB,RightB])
line([th_p,th_p],[0,max(N)])
line([th_n,th_n],[0,max(N)])

% x=linspace(min(EDGES),max(EDGES),10000);
% figure; hold on;plot(x,fitresult.a1*exp(-((x-fitresult.b1)./fitresult.c1).^2),'Color','r','Linewidth',2)

id=[1,2,4,5];%this is for DMDLDLCB leaving BR out for this analysis

for i=1:size(F_dff_pool_all,1)
    for j=1:length(id)
        if Fit_Exp_all(i,id(j))>=th_p
            dyn_code(i,j)=2;
        elseif Fit_Exp_all(i,id(j))<=th_n
            dyn_code(i,j)=1;
        elseif Fit_Exp_all(i,id(j))>=th_n & Fit_Exp_all(i,id(j))<=th_p
             dyn_code(i,j)=0;
        else
            dyn_code(i,j)=3;
        end
    end

end

dyn_code_m=dyn_code;
for i=1:size(dyn_code_m,1)
%     if ismember(3,dyn_code_m(i,:))%if after removing non-responsive stimuli still have no value for exp set the whole thing to 3
%         dyn_code_m(i,:)=[3,3,3,3];
%     end
    temp=num2str(dyn_code_m(i,:));
    temp=temp(find(~isspace(temp)));
    S(i)=base2dec(temp,4);
end
%% make a bar plot of the number of neurons in each cluster
for i=2:length(Xu)
    NU(i-1)=sum(X==Xu(i));
    index=find(X==Xu(i));
    temp=dec2bin(Xu(i),4);
    for j=1:4
        if (j==1 || j==2 || j==4) && temp(j)=='1'   %take neurons to respond either to DL or to CB or DM  
            Dcode_base(j)=1;
        else
            Dcode_base(j)=NaN;
        end
    end
%     Dcode=repmat(Dcode_base,NU(i-1),1);
%     PI(i-1)=nansum(reshape(dyn_code(X==Xu(i),:).*Dcode,[1, NU(i-1)*4])==2)/nansum(reshape(dyn_code(X==Xu(i),:).*Dcode,[1, NU(i-1)*4])>-1);
    Inc(i-1)=sum(dyn_code(X==Xu(i),1)==2 & temp(1)=='1' | dyn_code(X==Xu(i),2)==2 & temp(2)=='1'| dyn_code(X==Xu(i),3)==2 & temp(3)=='1' | dyn_code(X==Xu(i),4)==2 & temp(4)=='1');
    
    if temp(2)=='1'
        PT(i-1)=mean(trial_ave_ttc(X==Xu(i),2));
    elseif temp(4)=='1'
        PT(i-1)=mean(trial_ave_ttc(X==Xu(i),5));
    else
        PT(i-1)=NaN;
    end

end
%at least 10 units per fish in at least three fish
for i=1:length(Xu)
    temp=fn(X==Xu(i));
    fnu{i}=unique(temp);
    for k=1:length(fnu{i})
        unitcnt{i}(k)=sum(temp==fnu{i}(k));
    end
    %number of fish with more than 10 units in each
    fnucnt(i)=length(fnu{i}(unitcnt{i}>10))
end

%calculate average motor correlation for each of the clusters:
for i=1:length(Xu)
    N_cc(i)=100*sum(CCStimMot(X==Xu(i),6)>0.5)/length(CCStimMot(X==Xu(i),6));
    figure
    pie([N_cc(i),100-N_cc(i)])
    title(strcat('Cluster ',num2str(Xu(i))))
end
    
XuR=Xu(2:end);%responsive units
[B,I]=sort(NU,'descend');
figure;
bar(NU(I));
set(gca,'xticklabel',XuR(I))
hold on

load('binocularity_index.mat')
figure
BI2=BI(1:15);
bar(BI2(I));
set(gca,'xticklabel',XuR(I))
hold on


figure
[B,I]=sort(Inc,'descend');

bar(Inc(I)./NU(I));
set(gca,'xticklabel',XuR(I))
hold on



save('Peaks_EXPFits_ClustsNew_v2.mat')



