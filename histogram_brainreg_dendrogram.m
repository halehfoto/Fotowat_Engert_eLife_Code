%this is a program that reads the information about brain regions and
%creates histograms. Created by Haleh Fotowat
close all
clearvars
close all

path=uigetdir;
cd(path)
RegData=csvread('DMDLBRDLCB_celltypes_registered_v2.csv');
[NUMv,TXTv,RAWv]=xlsread('volume_mece.xlsx');
[NUM,TXT,RAW]=xlsread('DMDLBRDLCB_celltypes_registered_v2_Hierarical_allnodes.xlsx');
load('DMDLBRDLCB_082628_091112_HA.mat')
load('LS_DS_dynindices.mat')
ISI=40;
per2p=1.06295;
stim_type_cnt=5;
repse=5;
fps=69;
f0=round(ISI/(2*per2p));
nframes=size(F_dff_pool_all,2)
framerate=1/per2p;

time=1/framerate:1/framerate:nframes/framerate;
EI=RegData(:,7);
%% Plot a tree structure of all the aread in each three brain sections
All_regions_temp=unique(TXT(2:end,8));
All_regions=All_regions_temp(2:end);%the first row are the un-assigned cells/units
for i=1:length(All_regions)
  Regs=split(All_regions{i},',');
  ID=strmatch(Regs(end),TXTv(:,1)); %match the last node to the volume data
  All_Regs_vol(i)=sum(NUMv(ID-1)); %sometimes there is more than one volume data for a region (e.g. a region that spreads across two and the name contains a number which does not show up as character). combine the two volumes
end
k=1;
for i=2:length(TXT)
    if ~isempty(TXT{i,8})
        %Regs=split(TXT{i,8},',');
            %for j=1:length(Regs)
                ID=strmatch(TXT{i,8},All_regions);%find the cells that belong to any of the All_regions
                if ~isempty(ID)
                    Areas(k)=ID(1);%in case there were two matches in All_regions, take the first one as label (this is the same as described above with cross-regional areas identified with a *1)
                    fn_Area_Type(k)=fn(i-1);
                    k=k+1;
                end
            %end
    end
end
%find how many fish had cells in each brain region
for i=1:length(All_regions)
    All_fcnt(i)=length(unique(fn_Area_Type(Areas==i)));
    All_N(i)=length((fn_Area_Type(Areas==i)));
    if All_fcnt(i)<3 | All_N(i)<10 %only cells that are present in at least 2 fish
        All_fcnt(i)=0;
    end

end
    
[N,EDGES]=histcounts(Areas,'BinLimits',[1,length(All_regions)]);

% figure;bar(1:1:length(All_regions),N)

%normalize to the number of fish and the volume of the region.
Nn=N./(All_fcnt.*All_Regs_vol);
        
% figure;barh(1:1:length(All_regions),Nn)
% yticks(1:1:length(All_regions))
% set(gca,'yticklabel',All_regions)
% xtickangle(90)
% title('Wholebrain, all types')

figure;
index=find(~isnan(Nn) &~isinf(Nn)& Nn>0)
Nn2=Nn(index);

[B,I]=sort(Nn2,'descend');

for i=1:length(I)
    labels_temp=split(All_regions(index(I(i))),',');
    label(i)=labels_temp(end);
end
figure;bar(1:1:length(index(I)),Nn2((I)))
xticks(1:1:length(index(I)))
set(gca,'xticklabel',label)
xtickangle(90);


Type=[5,12];
for m=1:length(Type)
    k=1;
    for i=2:length(TXT)
        if ~isempty(TXT{i,8})
            %Regs=split(TXT{i,8},',');
            %for j=1:length(Regs)
            ID=strmatch(TXT{i,8},All_regions);
            if ~isempty(ID) & RegData(i-1,5)==Type(m)
                AreasR(k)=ID(1);
                fn_Area_TypeR(k)=fn(i-1);
                k=k+1;
            end
            %end
        end
    end
    %find how many fish had cells in each brain region
    for i=1:length(All_regions)
        All_fcntR(i)=length(unique(fn_Area_TypeR(AreasR==i)));
        N_total(i)=length((fn_Area_TypeR(AreasR==i)));
        if All_fcntR(i)<3 || N_total(i)<10 %only cells that are present in at least 3 fish and at least 10 cells
            All_fcntR(i)=0;
        end
        
    end
    [N,EDGES]=histcounts(AreasR,'BinLimits',[1,length(All_regions)]);
    %N(N<10)=0;%take only categories that had at least 10 cells
        
    %normalize to the number of fish and the volume of the region.
    
    NnR{m}=N./(All_fcntR.*All_Regs_vol)    
    clear AreasR
    clear fn_Area_TypeR
    clear N
    clear All_fcntR

end

NnA=[NnR{1};NnR{2}];
index=find((~isnan(NnR{1})&~isinf(NnR{1})&NnR{1}>0)|(~isnan(NnR{2})&~isinf(NnR{2})&NnR{2}>0));
N1=NnR{1}(index);
N1(isnan(N1))=0;
N1(isinf(N1))=0;

N2=NnR{2}(index);
N2(isnan(N2))=0;
N2(isinf(N2))=0;

[B,I]=sort(N1+N2,'descend');

for i=1:length(I)
    labels_temp=split(All_regions(index(I(i))),',');
    label(i)=labels_temp(end);
end
%Figure 5 d
figure;bar(1:1:length(index(I)),[NnR{2}(index(I));NnR{1}(index(I))],'stacked')
xticks(1:1:length(index(I)))
set(gca,'xticklabel',label)
xtickangle(90);

% figure;bar(1:1:length(index(I)),NnR{1}(index(I)),'stacked')
% xticks(1:1:length(index(I)))
% set(gca,'xticklabel',label)
% xtickangle(90);

label_all=label;
index_all=index;
I_all=I;
%% plot distribution of different dynamics for each cluster
%for DS neurons
X=RegData(:,5);
S=RegData(:,6);
load('Peaks_EXPFits_ClustsNew_v2.mat','dyn_code','expfit_signed','CCStimMot_th','th_p','th_n')
%calculate the decimal number based on a three digit code DM,DL,CB 
% for i=1:length(dyn_code)
%     temp=num2str(dyn_code(i,[1,2]));
%     temp=temp(find(~isspace(temp)));
%     S2(i)=base2dec(temp,4);
% 
% end
% e=min(S2):1:max(S2);
% SS=S2(X==12);
% combs=0:1:max(S2);
% combs_b4=dec2base(combs,4,2);
% % for i=1:length(combs_b4)
% %     combs_b4(i,end)='x';
% % end
% [N,EDGES]=histcounts(S2(X==12),16);
% figure;hold on;
% subplot(3,1,1)
% [B,I]=sort(N,'descend');
% bar(combs,100*N(I)/sum(N));
% set(gca,'xtick',combs)
% set(gca,'xticklabel',{[]})
% set(gca,'xticklabel',combs_b4(I,:))
% xtickangle(90)
% title('DS-DynamicClust')
% ylabel('percentage of neurons')
% 
% YLIM=ylim(gca)
% [N,EDGES]=histcounts(S2(X==14),16);
% hold on;
% subplot(3,1,2)
% [B,I]=sort(N,'descend');
% bar(combs,N(I)/(sum(N)));
% set(gca,'xtick',combs)
% set(gca,'xticklabel',{[]})
% set(gca,'xticklabel',combs_b4(I,:))
% xtickangle(90)
% title('DS2-DynamicClust')
% 
% 
% 
% clear S2
% %for LS neurons
% for i=1:length(dyn_code)
%     temp=num2str(dyn_code(i,[2,4]));
%     temp=temp(find(~isspace(temp)));
%     S2(i)=base2dec(temp,4);
% 
% end
% e=min(S2):1:max(S2);
% SS=S2(X==5);
% combs=0:1:max(S2);
% combs_b4=dec2base(combs,4,2);
% % for i=1:length(combs_b4)
% %     combs_b4(i,end)='x';
% % end
% [N,EDGES]=histcounts(S2(X==5),16);
% subplot(3,1,3)
% [B,I]=sort(N,'descend');
% bar(combs,100*N(I)/sum(N));
% set(gca,'xtick',combs)
% set(gca,'xticklabel',{[]})
% set(gca,'xticklabel',combs_b4(I,:))
% xtickangle(90)
% title('LS-DynamicClust')
% ylabel('percentage of neurons')



%% plot the relative distribution of the different groups New Figure 7
%% 1- LS
AreasR=[]
Areas_all=[]
fn_Area_TypeR=[]
All_fcntR=[]
All_fcnt_all=[]
fn_Area_Type_all=[]
clear N
clear NnA
dc=[1,0,2];
k=1;
for i=2:length(TXT)
    if ~isempty(TXT{i,8})
        %Regs=split(TXT{i,8},',');
        %for j=1:length(Regs)
        ID=strmatch(TXT{i,8},All_regions);
        if ~isempty(ID) & RegData(i-1,5)==5
            AreasR_X5(k)=ID(1);
            fn_Area_TypeR_X5(k)=fn(i-1);
            k=k+1;
        end
        %end
    end
end
%find how many fish had cells in each brain region
for i=1:length(All_regions)
    All_fcntR_X5(i)=length(unique(fn_Area_TypeR_X5(AreasR_X5==i)));
    N_X5=length((fn_Area_TypeR_X5(AreasR_X5==i)));
    if All_fcntR_X5(i)<3  || N_X5<10 %only cells that are present in at least 3 fish and at least 10 cells
        All_fcntR_X5(i)=0;
    end
    
end


for m=1:length(dc)
    k=1;
    kk=1;
    for i=2:length(TXT)
        if ~isempty(TXT{i,8})
            ID=strmatch(TXT{i,8},All_regions);
            if ~isempty(ID) & RegData(i-1,5)==5 & dyn_code(i-1,2)==dc(m)
                    AreasR(k)=ID(1);
                    fn_Area_TypeR(k)=fn(i-1);
                    k=k+1;
            end
            %end
        end
    end
    %find how many fish had cells in each brain region
    for i=1:length(All_regions)
        All_fcntR(i)=length(unique(fn_Area_TypeR(AreasR==i)));
%         All_fcnt_all(i)=length(unique(fn_Area_Type_all(Areas_all==i)));

        if All_fcntR_X5(i)==0  %only cells that are present in at least 3 fish and at least 10 cells
            All_fcntR(i)=0;
        end
        
    end
    [N,EDGES]=histcounts(AreasR,'BinLimits',[1,length(All_regions)]);
    %N(N<10)=0;%take only categories that had at least 10 cells
        
    %normalize to the number of fish and the volume of the region.
    
    NnR{m}=N./(All_fcntR.*All_Regs_vol)    
    clear AreasR
    clear fn_Area_TypeR

    clear N
    clear All_fcntR
    clear All_fcnt_all
        
end

NnA=[NnR{1};NnR{2};NnR{3}];
index=find((~isnan(NnR{1})&~isinf(NnR{1})&NnR{1}>0)|(~isnan(NnR{2})&~isinf(NnR{2})&NnR{2}>0)|(~isnan(NnR{3})&~isinf(NnR{3})&NnR{3}>0));
N1=NnR{1}(index);
N1(isnan(N1))=0;
N1(isinf(N1))=0;

N2=NnR{2}(index);
N2(isnan(N2))=0;
N2(isinf(N2))=0;

N3=NnR{3}(index);
N3(isnan(N3))=0;
N3(isinf(N3))=0;


[B,I]=sort(N1+N2+N3,'descend');
clear label
for i=1:length(I)
    labels_temp=split(All_regions(index(I(i))),',');
    label(i)=labels_temp(end);
end
% figure;H=bar(1:1:length(index(I)),[NnR{1}(index(I));NnR{2}(index(I));NnR{3}(index(I))],'stacked')
% 
% colors=[0 0 1; 0 1 1; 1 0 0]
% for i=1:3
%     H(i).FaceColor='flat';
%     H(i).CData= colors(i,:);
% end
% xticks(1:1:length(index(I)))
% set(gca,'xticklabel',label)
% xtickangle(90);
for i=1:length(index(I))
    labels_temp=split(All_regions(index(I(i))),',');
    label(i)=labels_temp(end);
end

figure;H=bar(1:1:length(index_all(I_all)),[NnR{1}(index_all(I_all));NnR{2}(index_all(I_all));NnR{3}(index_all(I_all))],'stacked')
colors=[0 0 1; 0 1 1; 1 0 0]
for i=1:3
    H(i).FaceColor='flat';
    H(i).CData= colors(i,:);
end

xticks(1:1:length(I_all))

set(gca,'xticklabel',label_all)
xtickangle(90)

%% 2- DS
% clear All_fcnt_all
% clear fn_Area_TypeR
% clear fn_Area_Type_all
% clear AreasR
% clear Areas_all
% clear All_fcntR
% clear All_fcnt_all
% 
k=1;
for i=2:length(TXT)
    if ~isempty(TXT{i,8})
        %Regs=split(TXT{i,8},',');
        %for j=1:length(Regs)
        ID=strmatch(TXT{i,8},All_regions);
        if ~isempty(ID) & RegData(i-1,5)==12
            AreasR_X12(k)=ID(1);
            fn_Area_TypeR_X12(k)=fn(i-1);
            k=k+1;
        end
        %end
    end
end
%find how many fish had cells in each brain region
for i=1:length(All_regions)
    All_fcntR_X12(i)=length(unique(fn_Area_TypeR_X12(AreasR_X12==i)));
    N_X12(i)=length((fn_Area_TypeR_X12(AreasR_X12==i)));
    if All_fcntR_X12(i)<3  || N_X12(i)<10 %only cells that are present in at least 3 fish and at least 10 cells
        All_fcntR_X12(i)=0;
    end
    
end

clear AreasR
clear fn_Area_TypeR
clear N
clear All_fcntR


dc=[1,0,2];
for m=1:length(dc)
    k=1;
    kk=1;
    for i=2:length(TXT)
        if ~isempty(TXT{i,8})
            %Regs=split(TXT{i,8},',');
            %for j=1:length(Regs)
            ID=strmatch(TXT{i,8},All_regions);
            if ~isempty(ID) & RegData(i-1,5)==12 & dyn_code(i-1,1)==dc(m)
                    AreasR(k)=ID(1);
                    fn_Area_TypeR(k)=fn(i-1);
                    k=k+1;
            end
            %end
        end
    end
    %find how many fish had cells in each brain region
    for i=1:length(All_regions)
        All_fcntR(i)=length(unique(fn_Area_TypeR(AreasR==i)));
        if All_fcntR_X12(i)==0  %only cells that are present in at least 3 fish and at least 10 cells in the whole group pooling all dynamic types
            All_fcntR(i)=0;
        end
        
    end
    [N,EDGES]=histcounts(AreasR,'BinLimits',[1,length(All_regions)]);
    %N(N<10)=0;%take only categories that had at least 10 cells
        
    %normalize to the number of fish and the volume of the region.
    
    NnR{m}=N./(All_fcntR.*All_Regs_vol)    
    clear AreasR
    clear fn_Area_TypeR
    clear N
    clear All_fcntR

end

NnA=[NnR{1};NnR{2};NnR{3}];
index=find((~isnan(NnR{1})&~isinf(NnR{1})&NnR{1}>0)|(~isnan(NnR{2})&~isinf(NnR{2})&NnR{2}>0)|(~isnan(NnR{3})&~isinf(NnR{3})&NnR{3}>0));
N1=NnR{1}(index);
N1(isnan(N1))=0;
N1(isinf(N1))=0;

N2=NnR{2}(index);
N2(isnan(N2))=0;
N2(isinf(N2))=0;

N3=NnR{3}(index);
N3(isnan(N3))=0;
N3(isinf(N3))=0;


[B,I]=sort(N1+N2+N3,'descend');
clear label
for i=1:length(I)
    labels_temp=split(All_regions(index(I(i))),',');
    label(i)=labels_temp(end);
end
% figure;H=bar(1:1:length(index(I)),[NnR{1}(index(I));NnR{2}(index(I));NnR{3}(index(I))],'stacked')
% colors=[0 0 1; 0 1 1; 1 0 0]
% for i=1:3
%     H(i).FaceColor='flat';
%     H(i).CData= colors(i,:);
% end
% xticks(1:1:length(index(I)))
% set(gca,'xticklabel',label)
% xtickangle(90);

figure;H=bar(1:1:length(index_all(I_all)),[NnR{1}(index_all(I_all));NnR{2}(index_all(I_all));NnR{3}(index_all(I_all))],'stacked')
colors=[0 0 1; 0 1 1; 1 0 0]
for i=1:3
    H(i).FaceColor='flat';
    H(i).CData= colors(i,:);
end

xticks(1:1:length(I_all))

set(gca,'xticklabel',label_all)
xtickangle(90)


%% plot the DS_I
figure;subplot(2,1,1);hold on;%plot(time, (F_dff_pool_all(X==12 & dyn_code(:,1)==2,:)),'Color',[0.8,0.8,0.8]); hold on
plot(time, median(F_dff_pool_all(X==12 & dyn_code(:,1)==2,:)/max(median(F_dff_pool_all(X==12 & dyn_code(:,1)==2,:)))),'Color','b','LineWidth',2);
xlim(time([1,fps*5]))
LIM=ylim(gca);

DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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
ylim(LIM)
% axis tight

title(strcat('DSI, n=',num2str(sum(X==12 & dyn_code(:,1)==2))))
subplot(2,1,2);%plot(time, (F_dff_pool_all(X==5 & dyn_code(:,2)==2,:)),'Color',[0.8,0.8,0.8]); hold on
hold on;plot(time, median(F_dff_pool_all(X==12 & dyn_code(:,2)==2,:)/max(median(F_dff_pool_all(X==12 & dyn_code(:,2)==2,:)))),'Color','r','LineWidth',2);

LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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
xlim(time([fps*5+1,fps*10]))
title(strcat('LSI, n=',num2str(sum(X==5 & dyn_code(:,2)==2))))
%% plot DSD LSD
figure;subplot(2,1,1);plot(time, median(F_dff_pool_all(X==12 & dyn_code(:,1)==1,:))/max(median(F_dff_pool_all(X==12 & dyn_code(:,1)==1,:))),'Color','b','LineWidth',2); hold on
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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
xlim(time([1,fps*5]))
title(strcat('DSD, n=',num2str(sum(X==12 & dyn_code(:,1)==1))))
subplot(2,1,2);plot(time, median(F_dff_pool_all(X==5 & dyn_code(:,2)==1,:))/max(median(F_dff_pool_all(X==5 & dyn_code(:,2)==1,:))),'Color','r','LineWidth',2); hold on
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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
xlim(time([fps*5+1,fps*10]))
title(strcat('LSD, n=',num2str(sum(X==5 & dyn_code(:,2)==1))))
%% plot DSNC LSNC
figure;subplot(2,1,1); hold on
%plot(time, (F_dff_pool_all(X==12 & dyn_code(:,1)==0,:)),'Color',[0.8,0.8,0.8]); 
plot(time, median(F_dff_pool_all(X==12 & dyn_code(:,1)==0,:))/max(median(F_dff_pool_all(X==12 & dyn_code(:,1)==0,:))),'Color','r','LineWidth',2);
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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
xlim(time([1,fps*5]))
title(strcat('DSNC, n=',num2str(sum(X==12 & dyn_code(:,1)==0))))
subplot(2,1,2); hold on
%plot(time, (F_dff_pool_all(X==5 & dyn_code(:,2)==0,:)),'Color',[0.8,0.8,0.8]); 
plot(time, median(F_dff_pool_all(X==5 & dyn_code(:,2)==0,:))/max(median(F_dff_pool_all(X==5 & dyn_code(:,2)==0,:))),'Color','b','LineWidth',2);
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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
xlim(time([fps*5+1,fps*10]))
title(strcat('LSNC, n=',num2str(sum(X==5 & dyn_code(:,1)==0))))

%% plot no fit
figure;subplot(2,1,1);plot(time, median(F_dff_pool_all(X==12 & dyn_code(:,2)==3,:))); hold on
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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
xlim(time([fps*5+1,fps*10]))
title('DSNF')
subplot(2,1,2);plot(time, median(F_dff_pool_all(X==5 & dyn_code(:,2)==3,:))); hold on
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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
xlim(time([fps*5+1,fps*10]))
title('LSNF')
%% plot example traces
%1100 DS
Type='3233';
C=12;
figure;stdshade(F_dff_pool_all(X==C& S==base2dec(Type,4),:),0.2,'k',time);
hold on
title(strcat('n=',num2str(sum(X==C& S==base2dec(Type,4))),', Type=',Type))
hold on
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'k','r','c','r','b'};

for k=1:5
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

EI_all=RegData(:,7);

dsi=[base2dec('1123',4),base2dec('1233',4),base2dec('2023',4),base2dec('2133',4),base2dec('2233',4)];
dsd=[base2dec('0133',4),base2dec('1033',4),base2dec('1333',4),base2dec('3133',4)];
sum(EI_all(X==12& ismember(S,dsi))==-1)/sum(EI_all(X==12& ismember(S,dsi))==1)
sum(EI_all(X==12& ismember(S,dsd))==-1)/sum(EI_all(X==12& ismember(S,dsd))==1)


