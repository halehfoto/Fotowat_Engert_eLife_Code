%this is a program that generates trial averages and forms clusters based
%on that/Select folder TwoEyeAnalysis/ISI20/LXLXCXCX
close all
clearvars
path=uigetdir;
cd(path)

%read the pooled clustering data
load('LXLXCXCX_HA.mat')
clear IDX

load('Peaks_EXPFits_LXLXCXCX_v2.mat')

stim_type_cnt=4;
repse=10;
per2p=1.06295;
ISI=20;
f0=round(ISI/(2*per2p));
nframes=size(F_dff_pool_all,2)
framerate=1/per2p;
fps=50;
time=1/framerate:1/framerate:nframes/framerate;




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
ax1=imagesc(trial_ave_mat(I,1:(fps+1)))
hold on
line([f0,f0],[1,length(I)],'color','w')
line([(fps+1)-f0,(fps+1)-f0],[1,length(I)],'color','w')
%colorbar
caxis([-2,6])

subplot(1,6,2)
%DL
ax2=imagesc(trial_ave_mat(I,(fps+1):(fps+1)*2))
line([f0,f0],[1,length(I)],'color','w')
line([(fps+1)-f0,(fps+1)-f0],[1,length(I)],'color','w')
%colorbar
caxis([-2,6])


subplot(1,6,3)
%CB
ax3=imagesc(trial_ave_mat(I,(fps+1)*2:(fps+1)*3)) %start 10 frames before stim onset
hold on
line([f0,f0],[1,length(I)],'color','w')
line([(fps+1)-f0,(fps+1)-f0],[1,length(I)],'color','w')
% colorbar
caxis([-2,6])


subplot(1,6,4)
%CB
ax4=imagesc(trial_ave_mat(I,(fps+1)*3:end)) %start 10 frames before stim onset
hold on
line([f0,f0],[1,length(I)],'color','w')
line([(fps+1)-f0,(fps+1)-f0],[1,length(I)],'color','w')
%colorbar
caxis([-2,6])

subplot(1,6,5)
%CB
ax5=imagesc(CCStimMot_th(I,:)) %start 10 frames before stim onset
hold on

ax6=subplot(1,6,6)
%CB
imagesc(X(I)');colormap(ax6,1-jet) %start 10 frames before stim onset
hold on

%plot the time course
figure
index=1:100:length(X);
colormap1=1-jet;

for i=1:length(index)-1
    %figure
    subplot(length(index)-1,1,i)
    hold on
    plot(1:(fps+1),median(trial_ave_mat(I(index(i):index(i+1)-1),1:(fps+1)),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:))
    line([f0,f0],[-0.5,0.5],'color','k')
    line([(fps+1)-f0,(fps+1)-f0],[-0.5,0.5],'color','k')
    plot((fps+1)+1:(fps+1)*2,median(trial_ave_mat(I(index(i):index(i+1)-1),(fps+1)+1:(fps+1)*2),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:))
    line([(fps+1)+f0,(fps+1)+f0],[-0.5,0.5],'color','k')
    line([(fps+1)*2-f0,(fps+1)*2-f0],[-0.5,0.5],'color','k')

    plot((fps+1)*2+1:(fps+1)*3,median(trial_ave_mat(I(index(i):index(i+1)-1),(fps+1)*2+1:(fps+1)*3),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:))
    line([(fps+1)*2+f0,(fps+1)*2+f0],[-0.5,0.5],'color','k')
    line([(fps+1)*3-f0,(fps+1)*3-f0],[-0.5,0.5],'color','k')
    
    plot((fps+1)*3+1:(fps+1)*4,median(trial_ave_mat(I(index(i):index(i+1)-1),(fps+1)*3+1:end),1),'Color',colormap1(17*(round(mean(X(I(index(i):index(i+1)-1)))))+1,:))
    line([(fps+1)*3+f0,(fps+1)*3+f0],[-0.5,0.5],'color','k')
    line([(fps+1)*4-f0,(fps+1)*4-f0],[-0.5,0.5],'color','k')
    %axis tight
    xlim([1,(fps+1)*4])
    ylim([-1,3])
    axis off
    
    
end

%define various response categories
Xu=unique(X);
for i=1:length(Xu)
    XuL(i)=sum(X==Xu(i));
    fnX(i)=length(unique(fn(X==Xu(i))));
end
%take only units that were present in at least 3 fish

%take only categories that contained more than 100 neurons
Index=find(XuL>50);
for i=2:length(Index)
    figure;
    
    title(strcat('Clust=',num2str(Xu(Index(i))),', n=',num2str(sum(X==Xu(Index(i))))))
    hold on
    stdshade(F_dff_pool_all(X==Xu(Index(i)),:),0.2,'k',time)
    hold on
    LIM=ylim(gca);
    DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
    kkk=1;
    color_bar={'r','m','b','c'};
    
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
end
colormap1=1-jet;
for i=2:length(Index)
    figure;
    title(strcat('Clust=',num2str(Xu(Index(i))),', n=',num2str(sum(X==Xu(Index(i))))))
    color=colormap1(Xu(Index(i))*17+1,:);
    hold on
    plot(time(f0:f0+fps),median(trial_ave_mat(X==Xu(Index(i)),1:(fps+1)),1),'Color',colormap1(17*Xu(Index(i))+1,:))
    hold on
    plot(time(f0+(fps+1):f0+(fps+1)*2),median(trial_ave_mat(X==Xu(Index(i)),fps+1:(fps+1)*2),1),'Color',colormap1(17*Xu(Index(i))+1,:))
    plot(time(f0+(fps+1)*2+1:f0+(fps+1)*3),median(trial_ave_mat(X==Xu(Index(i)),(fps+1)*2+1:(fps+1)*3),1),'Color',colormap1(17*Xu(Index(i))+1,:))
    plot(time(f0+(fps+1)*3+1:f0+(fps+1)*4),median(trial_ave_mat(X==Xu(Index(i)),(fps+1)*3+1:end),1),'Color',colormap1(17*Xu(Index(i))+1,:))
    axis tight
    LIM=ylim(gca);
    DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
    kkk=1;
    color_bar={'k','r','c','r','b'};
    
    for k=1:4
        %plot repse areas at first color
        for kk=1:1
          
            extra_time=5*per2p;
            area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
            area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
            
            kkk=kkk+1;
        end
        
    end
end
PeakInd_relon=NaN(size(PeakInd));
%compare response timing for the trial ave across different 
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

Peak_times=[trial_ave_ttc(X==1,2);trial_ave_ttc(X==2,2);trial_ave_ttc(X==8,2);trial_ave_ttc(X==4,2);trial_ave_ttc(X==12,2)];
groups=[ones(size(trial_ave_ttc(X==1,2)));2*ones(size(trial_ave_ttc(X==2,2)));3*ones(size(trial_ave_ttc(X==8,2)));4*ones(size(trial_ave_ttc(X==4,2)));5*ones(size(trial_ave_ttc(X==12,2)))];
[p,a,STATS]=kruskalwallis(Peak_times,groups);
multcompare(STATS)


IDX=X;
save('Peaks_EXPFits_ClustsNew_LXLXCXCX_v2.mat')

