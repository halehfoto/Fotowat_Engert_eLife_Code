%This is a program that plots different response categories in the two eye
%experiment. data is in TwoEyeAnalysis/ISI20/LXLXCXCX
close all
clearvars
path=uigetdir;
cd(path)
load('Peaks_EXPFits_ClustsNew_LXLXCXCX_v2.mat')

%% generate bar plot for all dimming responsive neurons

eye1=[0,1,2,3];
k=1;
for i=1:length(eye1)
    ratio_dim_all(k)=sum(IDX'==8 & dyn_code(:,1)==eye1(i) )/sum(IDX==12|IDX==8);
    k=k+1;
    figure;
    hold on
    index=find(IDX'==8 & dyn_code(:,1)==eye1(i));
    if length(index)>1
        
        stdshade(F_dff_pool_all(index,:),0.2,'k', time)
    end
    axis tight
    title(num2str(eye1(i)));
    %ylim([-1,3.5])
    LIM=ylim(gca);
    DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
    kkk=1;
    color_bar={'r','r','m','m'};
    for j=1:4
        %plot repse areas at first color
        for kk=1:10
            
            extra_time=5*per2p;
            area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{j},'FaceAlpha',0.1)
            area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{j},'FaceAlpha',0.1)
            
            kkk=kkk+1;
        end
        
    end
    
end
figure;
hold on
eye1=[1,2];

for i=1:length(eye1)
%     ratio_dim_all(k)=sum(IDX'==8 & dyn_code(:,1)==eye1(i) )/sum(IDX==12|IDX==8);
%     k=k+1;
    hold on
    index=find(IDX'==4 & dyn_code(:,2)==eye1(i));
    subplot(2,1,i)
    hold on
    if length(index)>1
        
        stdshade(F_dff_pool_all(index,1:20*fps),0.2,'c', time(1:20*fps))
    end
    axis tight
    title(num2str(eye1(i)));
    %ylim([-1,3.5])
    LIM=ylim(gca);
    DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
    kkk=1;
    color_bar={'r','m'};
    for j=1:2
        %plot repse areas at first color
        for kk=1:10
            
            extra_time=5*per2p;
            area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{j},'FaceAlpha',0.1)
            area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{j},'FaceAlpha',0.1)
            
            kkk=kkk+1;
        end
        
    end
    text(200,2,num2str(sum(IDX'==4 & dyn_code(:,2)==eye1(i))))
end

eye1=[0,1,2,3];
eye2=[0,1,2,3];
k=1;
for i=1:length(eye1)
    for j=1:length(eye2)
        ratio_bindim_all(k)=sum(IDX'==12 & dyn_code(:,1)==eye1(i) & dyn_code(:,2)==eye2(j))/sum(IDX==12|IDX==8);
        k=k+1;
        figure;
        hold on
        index=find(IDX'==12 & dyn_code(:,1)==eye1(i)&dyn_code(:,2)==eye2(j));
        if length(index)>1
            
            stdshade(F_dff_pool_all(index,:),0.2,'k', time)
        end
        axis tight
        title(strcat(num2str(eye1(i)),num2str(eye2(j))));
        
        %ylim([-1.5,3.5])
        LIM=ylim(gca);
        DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
        kkk=1;
        color_bar={'r','r','m','m'};
        for jj=1:4
            %plot repse areas at first color
            for kk=1:10
                
                extra_time=5*per2p;
                area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
                area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
                
                kkk=kkk+1;
            end
            
        end
        
        
    end
end
%% plot histogram of relative abundance
groups={'0';'1';'2';'3';'00';'01';'02';'03';'10';'11';'12';'13';'20';'21';'22';'23';'30';'31';'32';'33'};
dimall=[ratio_dim_all,ratio_bindim_all];
[I,B]=sort(dimall,'descend');
figure;bar(1:1:20,dimall(B))
xticks(1:20)
set(gca,'xticklabel',groups(B))
%% plot distribution of different clusters for all dimming responsive neuron response to the right eye
eye1=[1,2];
figure;
hold on
hold on

for i=1:length(eye1)
    subplot(2,1,i)
    hold on
    ratio_dim_R(i)=sum(IDX'==8 & dyn_code(:,1)==eye1(i) )/sum((IDX'==12|IDX'==8)& (dyn_code(:,1)==1|dyn_code(:,1)==2));
    ratio_bindim_R(i)=sum(IDX'==12 & dyn_code(:,1)==eye1(i) )/sum((IDX'==12|IDX'==8) & (dyn_code(:,1)==1|dyn_code(:,1)==2));
    if sum(IDX'==8 & dyn_code(:,1)==eye1(i))>1
        stdshade(F_dff_pool_all(IDX'==8 & dyn_code(:,1)==eye1(i),1:20*fps),0.2,'c',time(1:20*fps))
        LIM=ylim(gca);
        DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
        kkk=1;
        color_bar={'r','m'};
        for j=1:2
            %plot repse areas at first color
            for kk=1:10
                
                extra_time=5*per2p;
                area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{j},'FaceAlpha',0.1)
                area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{j},'FaceAlpha',0.1)
                
                kkk=kkk+1;
            end
            
        end
        text(500,2,num2str(sum(IDX'==8 & dyn_code(:,1)==eye1(i))))
        axis tight
    end
    
end
groups={'1X';'2X';' ';'1';'2'};

dim_all_sorted=[ratio_bindim_R(1),ratio_bindim_R(2),0,ratio_dim_R(1),ratio_dim_R(2)]
figure;bar(1:1:5,dim_all_sorted)
xticks(1:5)
set(gca,'xticklabel',groups)
%% do the same as above but pool the left and right side monocular dimming responsive
%% plot distribution of different clusters for all dimming responsive neuron response to the right eye
eye1=[1,2];
figure;
hold on
hold on

for i=1:length(eye1)
    subplot(2,1,i)
    hold on
    ratio_dim_both(i)=sum((IDX'==8 & dyn_code(:,1)==eye1(i))|(IDX'==4 & dyn_code(:,2)==eye1(i)) )/sum((IDX'==12 & (dyn_code(:,1)==1|dyn_code(:,1)==2))|(IDX'==8 & (dyn_code(:,1)==1|dyn_code(:,1)==2))|(IDX'==4 & (dyn_code(:,2)==1|dyn_code(:,2)==2)));
    ratio_bindim_R(i)=sum(IDX'==12 & dyn_code(:,1)==eye1(i) )/sum((IDX'==12 & (dyn_code(:,1)==1|dyn_code(:,1)==2))|(IDX'==8 & (dyn_code(:,1)==1|dyn_code(:,1)==2))|(IDX'==4 & (dyn_code(:,2)==1|dyn_code(:,2)==2)));
    if sum(IDX'==8 & dyn_code(:,1)==eye1(i))>1
        stdshade(F_dff_pool_all(IDX'==8 & dyn_code(:,1)==eye1(i),1:20*fps),0.2,'c',time(1:20*fps))
        LIM=ylim(gca);
        DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
        kkk=1;
        color_bar={'r','m'};
        for j=1:2
            %plot repse areas at first color
            for kk=1:10
                
                extra_time=5*per2p;
                area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{j},'FaceAlpha',0.1)
                area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{j},'FaceAlpha',0.1)
                
                kkk=kkk+1;
            end
            
        end
        text(500,2,num2str(sum(IDX'==8 & dyn_code(:,1)==eye1(i))))
        axis tight
    end
    
end
groups={'1';'2';' ';'1X';'2X'};

dim_all_sorted=[ratio_dim_both(1),ratio_dim_both(2),0,ratio_bindim_R(1),ratio_bindim_R(2)]
figure;bar(1:1:5,dim_all_sorted)
xticks(1:5)
set(gca,'xticklabel',groups)
%% plot the last trial first side to first trial second side trend for each of the three groups

index1=find(IDX'==12& dyn_code(:,1)==1);
%number of resetting and memory cells
N1r=0;
indexr=[];
N1m=0;
indexm=[];
figure;
hold on
for i=1:length(index1)
    if (mean(PeakAmps(index1(i),2,1:5))<=mean(PeakAmps(index1(i),1,6:10)))% & (dyn_code(index1(i),2)==1 | dyn_code(index1(i),2)==0)
        color=[0 0.8 1];
        N1m=N1m+1;
        indexm=[indexm,index1(i)];
    else
        if dyn_code(index1(i),2)==1
            color=[0.8 0.8 0.8];
            N1r=N1r+1;
            indexr=[indexr,index1(i)];
        end
    end
    plot([PeakAmps(index1(i),1,10);PeakAmps(index1(i),2,1)],'Color',color,'LineWidth',0.1)
end
plot([median(PeakAmps(indexm,1,10));median(PeakAmps(indexm,2,1))],'b','LineWidth',2)
plot([median(PeakAmps(indexr,1,10));median(PeakAmps(indexr,2,1))],'k','LineWidth',2)
figure;
subplot(2,1,1)
hold on
stdshade(F_dff_pool_all(indexm,1:20*fps),0.2,'b', time(1:20*fps))
LIM=ylim(gca);

kkk=1;
color_bar={'r','m'};
for jj=1:2
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight
text(1000,3,num2str(length(indexm)))
subplot(2,1,2)
hold on
stdshade(F_dff_pool_all(indexr,1:20*fps),0.2,'b', time(1:20*fps))
LIM=ylim(gca);

kkk=1;
color_bar={'r','m'};
for jj=1:2
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight
text(1000,3,num2str(length(indexr)))
%%


    
        
        
        
index2=find(IDX'==12& dyn_code(:,1)==2);
%number of resetting and memory cells
N2r=0;
index2r=[];
N2m=0;
index2m=[];
figure;
hold on
for i=1:length(index2)
    if mean(PeakAmps(index2(i),2,1:5))>=mean(PeakAmps(index2(i),1,6:10)) 
        color=[1 0.6 0];
        
        N2m=N2m+1;
        index2m=[index2m,index2(i)];
    else
        if dyn_code(index2(i),2)==2
            color=[0.8 0.8 0.8];
            N2r=N2r+1;
            index2r=[index2r,index2(i)];
        end
    end
    plot([PeakAmps(index2(i),1,10);PeakAmps(index2(i),2,1)],'Color',color,'LineWidth',0.1)
end
plot([median(PeakAmps(index2m,1,10));median(PeakAmps(index2m,2,1))],'r','LineWidth',2)
plot([median(PeakAmps(index2r,1,10));median(PeakAmps(index2r,2,1))],'k','LineWidth',2)

figure;
subplot(2,1,1)
hold on
stdshade(F_dff_pool_all(index2m,1:20*fps),0.2,'b', time(1:20*fps))
LIM=ylim(gca);

kkk=1;
color_bar={'r','m'};
for jj=1:2
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight
text(1000,2,num2str(length(index2m)))
subplot(2,1,2)
hold on
stdshade(F_dff_pool_all(index2r,1:20*fps),0.2,'b', time(1:20*fps))
LIM=ylim(gca);

kkk=1;
color_bar={'r','m'};
for jj=1:2
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{jj},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight
text(1000,2,num2str(length(index2r)))