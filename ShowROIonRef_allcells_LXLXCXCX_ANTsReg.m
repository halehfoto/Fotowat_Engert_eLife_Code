%this is a program to plot the ROIs on the reference brain
%read the reference brain
close all
clearvars
path=uigetdir;
cd(path)

filenames=dir('2*AllBrains_procdata_ha.mat');
for i=1:length(filenames)
    load(filenames(i).name,'Exp_seq')
    Eye{i}=Exp_seq(2);
end
[X,meta]=nrrdread('Elavl3-H2BRFP.nrrd');


%read the pooled clustering data
load('LXLXCXCX_HA.mat')
stim_type_cnt=4;
repse=10;
per2p=1.06295;
ISI=20;
f0=round(ISI/(2*per2p));
nframes=size(F_dff_pool_all,2)
framerate=1/per2p;
fps=50;
time=1/framerate:1/framerate:nframes/framerate;

RegData=csvread('LXLXCXCX_celltypes_registered_v2.csv');

fish_n=RegData(:,1)
% [NUM,TXT,RAW]=xlsread('DMDLBRDLCB_celltypes_registered_distribution.xlsx');
% colors=hsv(ccount);
%read the transformed ROIs

Y=std(single(X(:,:,:)),[],3);
Z=std(single(X(:,:,:)),[],1);
temp=(squeeze(Z))';
Zint=interp2(repmat(1:1:size(temp,1),size(temp,2),1),(repmat((1:1:size(temp,2))',1,size(temp,1))),temp',repmat(0.78/2:0.78/2:size(temp,1),size(temp,2),1),(repmat((1:1:size(temp,2))',1,length(0.78/2:0.78/2:size(temp,1)))));
Clusts=unique(RegData(:,5));
IDX=RegData(:,5);
%%
ROI_pix=(RegData(:,2:4))./([0.798,0.798,0.798]);
Zplanes=unique(round(ROI_pix(:,3)));   
colors=jet(length(Clusts));


%%
%claculate the binocularity index
%calculate binocularity index for dimming sensitive neurons
%i.e. clusters 8,4 and 12.
%first find indices for all dimming responsive neurons that respond to the
%right eye stimulation
index_dim_R=[];
index_dim_L=[];
for i=1:length(IDX)
    if (IDX(i)==8 & strmatch(Eye{fn(i)},'R'))| (IDX(i)==4 & strmatch(Eye{fn(i)},'L'))
        index_dim_R=[index_dim_R,i];
    elseif (IDX(i)==8 & strmatch(Eye{fn(i)},'L'))| (IDX(i)==4 & strmatch(Eye{fn(i)},'R'))
        index_dim_L=[index_dim_L,i];
    end
end
index_dim_R1=[];
for i=1:length(IDX)
    if (IDX(i)==8 & strmatch(Eye{fn(i)},'R'))
        index_dim_R1=[index_dim_R1,i];
    end
end
index_dim_R2=[];
for i=1:length(IDX)
    if (IDX(i)==4 & strmatch(Eye{fn(i)},'L'))
        index_dim_R2=[index_dim_R2,i];
    end
end


index_dim_L1=[];
for i=1:length(IDX)
    if (IDX(i)==8 & strmatch(Eye{fn(i)},'L'))
        index_dim_L1=[index_dim_L1,i];
    end
end
index_dim_L2=[];
for i=1:length(IDX)
    if (IDX(i)==4 & strmatch(Eye{fn(i)},'R'))
        index_dim_L2=[index_dim_L2,i];
    end
end
load('Peaks_EXPFits_ClustsNew_LXLXCXCX_v2.mat')

%% plot the right eye responses and contra side

index_dimR1_I=index_dim_R1(ROI_pix(index_dim_R1,2)<310.5);
index_dimR1_C=index_dim_R1(ROI_pix(index_dim_R1,2)>310.5);

index_dimR2_I=index_dim_R2(ROI_pix(index_dim_R2,2)<310.5);
index_dimR2_C=index_dim_R2(ROI_pix(index_dim_R2,2)>310.5);

%calculate %increasing neurons in each of the four categories
R1_I_inc=sum(dyn_code(index_dimR1_I,1)==2)/length(index_dimR1_I)
R2_I_inc=sum(dyn_code(index_dimR2_I,2)==2)/length(index_dimR2_I)



R1_C_inc=sum(dyn_code(index_dimR1_C,1)==2)/length(index_dimR1_C)
R2_C_inc=sum(dyn_code(index_dimR2_C,2)==2)/length(index_dimR2_C)


figure;plot(time,mean(F_dff_pool_all(index_dim_R1,:)));hold on
%ipsilateral
figure;plot(time, mean(F_dff_pool_all(index_dim_R1(ROI_pix(index_dim_R1,2)<310.5),:)),'r');hold on
%contralateral
plot(time, mean(F_dff_pool_all(index_dim_R1(ROI_pix(index_dim_R1,2)>310.5),:)),'g')
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};

BI_dim_R=(sum(ROI_pix(index_dim_R,2)<310.5)-sum(ROI_pix(index_dim_R,2)>310.5))/(sum(ROI_pix(index_dim_R,2)<310.5)+sum(ROI_pix(index_dim_R,2)>310.5))

for k=1:4
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight

figure;plot(time,mean(F_dff_pool_all(index_dim_R2,:)));hold on
%ipsilateral
figure;plot(time, mean(F_dff_pool_all(index_dim_R2(ROI_pix(index_dim_R2,2)<310.5),:)),'r');hold on
%contralateral
plot(time, mean(F_dff_pool_all(index_dim_R2(ROI_pix(index_dim_R2,2)>310.5),:)),'g')
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};


for k=1:4
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight


%% do the same for left eye
figure;plot(time,mean(F_dff_pool_all(index_dim_L1,:)));hold on
%ipsilateral
index_dimL1_I=index_dim_L1(ROI_pix(index_dim_L1,2)>310.5);
index_dimL1_C=index_dim_L1(ROI_pix(index_dim_L1,2)<310.5);

index_dimL2_I=index_dim_L2(ROI_pix(index_dim_L2,2)>310.5);
index_dimL2_C=index_dim_L2(ROI_pix(index_dim_L2,2)<310.5);

%calculate %increasing neurons in each of the four categories
L1_I_inc=sum(dyn_code(index_dimL1_I,1)==2)/length(index_dimL1_I)
L2_I_inc=sum(dyn_code(index_dimL2_I,2)==2)/length(index_dimL2_I)



L1_C_inc=sum(dyn_code(index_dimL1_C,1)==2)/length(index_dimL1_C)
L2_C_inc=sum(dyn_code(index_dimL2_C,2)==2)/length(index_dimL2_C)


figure;plot(time, mean(F_dff_pool_all(index_dim_L1(ROI_pix(index_dim_L1,2)>310.5),:)),'r');hold on
%contralateral
plot(time, mean(F_dff_pool_all(index_dim_L1(ROI_pix(index_dim_L1,2)<310.5),:)),'g')
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};


for k=1:4
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight

figure;plot(time,mean(F_dff_pool_all(index_dim_L2,:)));hold on
%ipsilateral
figure;plot(time, mean(F_dff_pool_all(index_dim_L2(ROI_pix(index_dim_L2,2)>310.5),:)),'r');hold on
%contralateral
plot(time, mean(F_dff_pool_all(index_dim_L2(ROI_pix(index_dim_L2,2)<310.5),:)),'g')
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};


for k=1:4
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight

%%

figure;imagesc(Y);colormap gray;hold on;axis equal
plot(ROI_pix(index_dim_R,1),ROI_pix(index_dim_R,2),'Color','c','LineStyle','none','Marker','.','MarkerSize',10)
plot(ROI_pix(index_dim_L,1),ROI_pix(index_dim_L,2),'Color','m','LineStyle','none','Marker','.','MarkerSize',10)


index_loom_R=[];
index_loom_L=[];

for i=1:length(IDX)
    if (IDX(i)==10 & strmatch(Eye{fn(i)},'R'))| (IDX(i)==5 & strmatch(Eye{fn(i)},'L'))
        index_loom_R=[index_loom_R,i];
    elseif (IDX(i)==10 & strmatch(Eye{fn(i)},'L'))| (IDX(i)==5 & strmatch(Eye{fn(i)},'R'))
        index_loom_L=[index_loom_L,i];
    end
end
        
figure;imagesc(Y);colormap gray;hold on;axis equal
plot(ROI_pix(index_loom_R,1),ROI_pix(index_loom_R,2),'Color','r','LineStyle','none','Marker','.','MarkerSize',10)
plot(ROI_pix(index_loom_L,1),ROI_pix(index_loom_L,2),'Color',[1,0.5,0.1],'LineStyle','none','Marker','.','MarkerSize',10)


index_bindim_R=[];%from the set of fish where the right eye was stimulated first
index_bindim_L=[];%from the set of fish where the left eye was stimulated first
for i=1:length(IDX)
    if (IDX(i)==12 & strmatch(Eye{fn(i)},'R'))
        index_bindim_R=[index_bindim_R,i];
    elseif (IDX(i)==12 & strmatch(Eye{fn(i)},'L'))
        index_bindim_L=[index_bindim_L,i];
    end
end
index_bindim_all=[index_bindim_R,index_bindim_L];
        
figure;imagesc(Y);colormap gray;hold on;axis equal
plot(ROI_pix(index_bindim_all,1),ROI_pix(index_bindim_all,2),'Color','b','LineStyle','none','Marker','.','MarkerSize',10)
% plot(ROI_pix(index_bindim_L,1),ROI_pix(index_bindim_L,2),'Color','c','LineStyle','none','Marker','.','MarkerSize',10)
%% separate by responsiveness to left or right eye
idx_rightr1=[];
idx_leftr1=[];
idx_rightl1=[];
idx_leftl1=[];

% for i=1:length(IDX)
%     if (IDX(i)==12 & strmatch(Eye{fn(i)},'R'))
%         if nanmean(PeakAmps(i,1,:))>=nanmean(PeakAmps(i,2,:))
%             idx_rightr1=[idx_rightr1,i];
%         else
%             idx_leftr1=[idx_leftr1,i];
%         end
%     elseif (IDX(i)==12 & strmatch(Eye{fn(i)},'L'))
%         if nanmean(PeakAmps(i,1,:))>=nanmean(PeakAmps(i,2,:))
%             idx_rightl1=[idx_rightl1,i];
%         else
%             idx_leftl1=[idx_leftl1,i];
%         end
%     end
% end
%compare only the first trial
for i=1:length(IDX)
    if (IDX(i)==12 & strmatch(Eye{fn(i)},'R'))
        if (PeakAmps(i,1,:))>=(PeakAmps(i,2,1))
            idx_rightr1=[idx_rightr1,i];
        else
            idx_leftr1=[idx_leftr1,i];
        end
    elseif (IDX(i)==12 & strmatch(Eye{fn(i)},'L'))
        if (PeakAmps(i,1,1))<(PeakAmps(i,2,1))
            idx_rightl1=[idx_rightl1,i];
        else
            idx_leftl1=[idx_leftl1,i];
        end
    end
end

figure;imagesc(Y);colormap gray;hold on;axis equal;
plot(ROI_pix(idx_rightr1,1),ROI_pix(idx_rightr1,2),'Color','b','LineStyle','none','Marker','.','MarkerSize',10)
plot(ROI_pix(idx_leftr1,1),ROI_pix(idx_leftr1,2),'Color','g','LineStyle','none','Marker','.','MarkerSize',10)
   
BI_bindim_R_r1=(sum(ROI_pix(idx_rightr1,2)<310.5)-sum(ROI_pix(idx_rightr1,2)>310.5))/(sum(ROI_pix(idx_rightr1,2)<310.5)+sum(ROI_pix(idx_rightr1,2)>310.5))
BI_bindim_L_r1=(sum(ROI_pix(idx_leftr1,2)>310.5)-sum(ROI_pix(idx_leftr1,2)<310.5))/(sum(ROI_pix(idx_leftr1,2)<310.5)+sum(ROI_pix(idx_leftr1,2)>310.5))

figure;imagesc(Y);colormap gray;hold on;axis equal
plot(ROI_pix(idx_rightl1,1),ROI_pix(idx_rightl1,2),'Color','b','LineStyle','none','Marker','.','MarkerSize',10)
plot(ROI_pix(idx_leftl1,1),ROI_pix(idx_leftl1,2),'Color','g','LineStyle','none','Marker','.','MarkerSize',10)

BI_bindim_R_l1=(sum(ROI_pix(idx_rightl1,2)<310.5)-sum(ROI_pix(idx_rightl1,2)>310.5))/(sum(ROI_pix(idx_rightl1,2)<310.5)+sum(ROI_pix(idx_rightl1,2)>310.5))
BI_bindim_L_l1=(sum(ROI_pix(idx_leftl1,2)>310.5)-sum(ROI_pix(idx_leftl1,2)<310.5))/(sum(ROI_pix(idx_leftl1,2)<310.5)+sum(ROI_pix(idx_leftl1,2)>310.5))

%calculate binocularity index
%for right eye stimulation
%% plot average responses for these clusters

BI_bindim_R=(sum(ROI_pix(index_bindim_R,2)<310.5)-sum(ROI_pix(index_bindim_R,2)>310.5))/(sum(ROI_pix(index_bindim_R,2)<310.5)+sum(ROI_pix(index_bindim_R,2)>310.5))
BI_bindim_L=(sum(ROI_pix(index_bindim_L,2)>310.5)-sum(ROI_pix(index_bindim_L,2)<310.5))/(sum(ROI_pix(index_bindim_L,2)<310.5)+sum(ROI_pix(index_bindim_L,2)>310.5))

%calculate the idex relative to the right eye
BI_bindim=(sum(ROI_pix(index_bindim_all,2)<310.5)-sum(ROI_pix(index_bindim_all,2)>310.5))/(sum(ROI_pix(index_bindim_all,2)<310.5)+sum(ROI_pix(index_bindim_all,2)>310.5))



%BI_bindim_R=sum(ROI_pix(index_bindim_all,2)>310.5)/sum(ROI_pix(index_bindim_all,2)<310.5)

BI_dim_R=(sum(ROI_pix(index_dim_R,2)<310.5)-sum(ROI_pix(index_dim_R,2)>310.5))/(sum(ROI_pix(index_dim_R,2)<310.5)+sum(ROI_pix(index_dim_R,2)>310.5))
BI_dim_L=(sum(ROI_pix(index_dim_L,2)>310.5)-sum(ROI_pix(index_dim_L,2)<310.5))/(sum(ROI_pix(index_dim_L,2)<310.5)+sum(ROI_pix(index_dim_L,2)>310.5))

figure;
BI_loom_R=(sum(ROI_pix(index_loom_R,2)<310.5)-sum(ROI_pix(index_loom_R,2)>310.5))/(sum(ROI_pix(index_loom_R,2)<310.5)+sum(ROI_pix(index_loom_R,2)>310.5))
BI_loom_L=(sum(ROI_pix(index_loom_L,2)>310.5)-sum(ROI_pix(index_loom_L,2)<310.5))/(sum(ROI_pix(index_loom_L,2)<310.5)+sum(ROI_pix(index_loom_L,2)>310.5))

%dimming responsive to the Right eye stimulation

figure
stdshade([trial_ave_mat(index_dim_R(IDX(index_dim_R)==8),1:(fps+1));trial_ave_mat(index_dim_R(IDX(index_dim_R)==4),(fps+1)+1:(fps+1)*2)],1,'b',time(1:(fps+1)))
hold on
stdshade([trial_ave_mat(index_dim_R(IDX(index_dim_R)==8),(fps+1)+1:(fps+1)*2);trial_ave_mat(index_dim_R(IDX(index_dim_R)==4),1:(fps+1))],1,'b',time((fps+1)+1:(fps+1)*2))


stdshade([trial_ave_mat(index_dim_R(IDX(index_dim_R)==8),(fps+1)*2+1:(fps+1)*3);trial_ave_mat(index_dim_R(IDX(index_dim_R)==4),(fps+1)*3+1:end)],1,'b',time((fps+1)*2+1:(fps+1)*3))
hold on
stdshade([trial_ave_mat(index_dim_R(IDX(index_dim_R)==8),(fps+1)*3+1:end);trial_ave_mat(index_dim_R(IDX(index_dim_R)==4),(fps+1)*2+1:(fps+1)*3)],1,'b',time((fps+1)*3+1:(fps+1)*4))
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};

for k=1:4
    %plot repse areas at first color
    for kk=1:10
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight

%dimming responsive to the left eye stimulation

figure
stdshade([trial_ave_mat(index_dim_L(IDX(index_dim_L)==4),1:(fps+1));trial_ave_mat(index_dim_L(IDX(index_dim_L)==8),(fps+1)+1:(fps+1)*2)],1,'b',time(1:(fps+1)))
hold on
stdshade([trial_ave_mat(index_dim_L(IDX(index_dim_L)==4),(fps+1)+1:(fps+1)*2);trial_ave_mat(index_dim_L(IDX(index_dim_L)==8),1:(fps+1))],1,'b',time((fps+1)+1:(fps+1)*2))


stdshade([trial_ave_mat(index_dim_L(IDX(index_dim_L)==4),(fps+1)*2+1:(fps+1)*3);trial_ave_mat(index_dim_L(IDX(index_dim_L)==8),(fps+1)*3+1:end)],1,'b',time((fps+1)*2+1:(fps+1)*3))
hold on
stdshade([trial_ave_mat(index_dim_L(IDX(index_dim_L)==4),(fps+1)*3+1:end);trial_ave_mat(index_dim_L(IDX(index_dim_L)==8),(fps+1)*2+1:(fps+1)*3)],1,'b',time((fps+1)*3+1:(fps+1)*4))
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};

for k=1:4
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight




%% Looming responsive to the Right eye stimulation
figure
stdshade([trial_ave_mat(index_loom_R(IDX(index_loom_R)==10),1:(fps+1));trial_ave_mat(index_loom_R(IDX(index_loom_R)==5),(fps+1)+1:(fps+1)*2)],1,'r',time(1:(fps+1)))
hold on
stdshade([trial_ave_mat(index_loom_R(IDX(index_loom_R)==10),(fps+1)+1:(fps+1)*2);trial_ave_mat(index_loom_R(IDX(index_loom_R)==5),1:(fps+1))],1,'r',time((fps+1)+1:(fps+1)*2))


stdshade([trial_ave_mat(index_loom_R(IDX(index_loom_R)==10),(fps+1)*2+1:(fps+1)*3);trial_ave_mat(index_loom_R(IDX(index_loom_R)==5),(fps+1)*3+1:end)],1,'r',time((fps+1)*2+1:(fps+1)*3))
hold on
stdshade([trial_ave_mat(index_loom_R(IDX(index_loom_R)==10),(fps+1)*3+1:end);trial_ave_mat(index_loom_R(IDX(index_loom_R)==5),(fps+1)*2+1:(fps+1)*3)],1,'r',time((fps+1)*3+1:(fps+1)*4))
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};

for k=1:4
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight


%% Looming responsive to the left eye stimulation

figure
stdshade([trial_ave_mat(index_loom_L(IDX(index_loom_L)==5),1:(fps+1));trial_ave_mat(index_loom_L(IDX(index_loom_L)==10),(fps+1)+1:(fps+1)*2)],1,'r',time(1:(fps+1)))
hold on
stdshade([trial_ave_mat(index_loom_L(IDX(index_loom_L)==5),(fps+1)+1:(fps+1)*2);trial_ave_mat(index_loom_L(IDX(index_loom_L)==10),1:(fps+1))],1,'r',time((fps+1)+1:(fps+1)*2))


stdshade([trial_ave_mat(index_loom_L(IDX(index_loom_L)==5),(fps+1)*2+1:(fps+1)*3);trial_ave_mat(index_loom_L(IDX(index_loom_L)==10),(fps+1)*3+1:end)],1,'r',time((fps+1)*2+1:(fps+1)*3))
hold on
stdshade([trial_ave_mat(index_loom_L(IDX(index_loom_L)==5),(fps+1)*3+1:end);trial_ave_mat(index_loom_L(IDX(index_loom_L)==10),(fps+1)*2+1:(fps+1)*3)],1,'r',time((fps+1)*3+1:(fps+1)*4))

LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};

for k=1:4
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight
%% binocular dimming stimulation all fish pooled
figure
stdshade(trial_ave_mat(index_bindim_all,1:(fps+1)),1,'b',time(1:(fps+1)))
hold on
stdshade(trial_ave_mat(index_bindim_all,(fps+1)+1:(fps+1)*2),1,'b',time((fps+1)+1:(fps+1)*2))


stdshade(trial_ave_mat(index_bindim_all,(fps+1)*2+1:(fps+1)*3),1,'b',time((fps+1)*2+1:(fps+1)*3))
hold on
stdshade(trial_ave_mat(index_bindim_all,(fps+1)*3+1:end),1,'b',time((fps+1)*3+1:(fps+1)*4))
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};

for k=1:4
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight

%% binocular dimming right eye stim first (different set of fish)
figure
stdshade(trial_ave_mat(index_bindim_R,1:(fps+1)),1,'b',time(1:(fps+1)))
hold on
stdshade(trial_ave_mat(index_bindim_R,(fps+1)+1:(fps+1)*2),1,'b',time((fps+1)+1:(fps+1)*2))


stdshade(trial_ave_mat(index_bindim_R,(fps+1)*2+1:(fps+1)*3),1,'b',time((fps+1)*2+1:(fps+1)*3))
hold on
stdshade(trial_ave_mat(index_bindim_R,(fps+1)*3+1:end),1,'b',time((fps+1)*3+1:(fps+1)*4))
LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};

for k=1:4
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight

%% binocular dimming left eye stim first (different set of fish)
figure
stdshade(trial_ave_mat(index_bindim_L,1:(fps+1)),1,'b',time(1:(fps+1)))
hold on
stdshade(trial_ave_mat(index_bindim_L,(fps+1)+1:(fps+1)*2),1,'b',time((fps+1)+1:(fps+1)*2))


stdshade(trial_ave_mat(index_bindim_L,(fps+1)*2+1:(fps+1)*3),1,'b',time((fps+1)*2+1:(fps+1)*3))
hold on
stdshade(trial_ave_mat(index_bindim_L,(fps+1)*3+1:end),1,'b',time((fps+1)*3+1:(fps+1)*4))

LIM=ylim(gca);
DT=time(fps-5)-ISI;%5 extra 2p frames at the end of each stimulus
kkk=1;
color_bar={'r','r','m','m'};

for k=1:4
    %plot repse areas at first color
    for kk=1:1
        
        extra_time=5*per2p;
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(2),LIM(2)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        area([kkk*ISI+(kkk-1)*DT+(kkk-1)*extra_time,kkk*(ISI+DT)+(kkk-1)*extra_time],[LIM(1),LIM(1)],'EdgeColor','None','FaceColor',color_bar{k},'FaceAlpha',0.1)
        
        kkk=kkk+1;
    end
    
end
axis tight








hold on
for i=1:length(index)
    if index(i)==12
        color='b'
    else
        color='r'
    end
%     subplot(length(index),1,i)
%     hold on
    stdshade(trial_ave_mat(X==index(i),1:(fps+1)),1,color,time(1:(fps+1)))
%     line([f0,f0],[-0.5,0.5],'color','k')
%     line([(fps+1)-f0,(fps+1)-f0],[-0.5,0.5],'color','k')
    
    stdshade(trial_ave_mat(X==index(i),(fps+1)+1:(fps+1)*2),1,color,time((fps+1)+1:(fps+1)*2))
%     line([(fps+1)+f0,(fps+1)+f0],[-0.5,0.5],'color','k')
%     line([(fps+1)*2-f0,(fps+1)*2-f0],[-0.5,0.5],'color','k')

    stdshade(trial_ave_mat(X==index(i),(fps+1)*2+1:(fps+1)*3),1,color,time((fps+1)*2+1:(fps+1)*3))
%     line([(fps+1)*2+f0,(fps+1)*2+f0],[-0.5,0.5],'color','k')
%     line([(fps+1)*3-f0,(fps+1)*3-f0],[-0.5,0.5],'color','k')
    
    stdshade(trial_ave_mat(X==index(i),(fps+1)*3+1:end),1,color,time((fps+1)*3+1:(fps+1)*4))
%     line([(fps+1)*3+f0,(fps+1)*3+f0],[-0.5,0.5],'color','k')
%     line([(fps+1)*4-f0,(fps+1)*4-f0],[-0.5,0.5],'color','k')
    %axis tight
    xlim([1,(fps+1)*4])
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

%plot LS neurons together (cluster 5 and 10)

figure(100)
imagesc(Y);
colormap gray
hold on
axis equal
figure(101)
%figure(2)
%imagesc((squeeze(Z))');
imagesc(Zint');
colormap gray
hold on
set(gca,'YDir','normal')
axis equal
axis tight
colormap1=1-jet;
for i=5
    
        index=find(IDX==i & flag==1);
        for k=1:length(index)
            if strmatch(Eye{fn(index(k))},'R')
                if i==10
                    cl='r';
                elseif i==5
                    cl='m';
                end
                %cl='r';
            else
                if i==10
                    cl='m';
                elseif i==5
                    cl='r';
                    
                end
                %cl='b';
            end
                
            figure(100)
            hold on
            title(strcat('Cluster',num2str(i)))
            %figure(1)
            plot(ROI_pix(index(k),1),ROI_pix(index(k),2),'Color',cl,'Marker','.','Markersize',6);
            figure(101)
            hold on
            title(strcat('Cluster',num2str(i)))

            %figure(2)
            plot(ROI_pix(index(k),1),ROI_pix(index(k),3),'Color',cl,'Marker','.','Markersize',6);
        end
    
    
end

figure(200)
imagesc(Y);
colormap gray
hold on
axis equal
figure(201)
%figure(2)
%imagesc((squeeze(Z))');
imagesc(Zint');
colormap gray
hold on
set(gca,'YDir','normal')
axis equal
axis tight
colormap1=1-jet;
for i=12
    
        index=find(IDX==i & flag==1);
        for k=1:length(index)
                
            figure(200)
            hold on
            title(strcat('Cluster',num2str(i)))
            %figure(1)
            plot(ROI_pix(index(k),1),ROI_pix(index(k),2),'Color','b','Marker','s','Markersize',6);
            figure(201)
            hold on
            title(strcat('Cluster',num2str(i)))

            %figure(2)
            plot(ROI_pix(index(k),1),ROI_pix(index(k),3),'Color','b','Marker','s','Markersize',6);
        end
    
    
end
%% generate bar plot for all dimming responsive neurons
eye1=[0,1,2,3];
eye2=[0,1,2,3];
k=1;
for i=1:length(eye1)
    for j=1:length(eye2)
        ratio_bindim_all(k)=sum(IDX'==12 & dyn_code(:,1)==eye1(i) & dyn_code(:,2)==eye2(j))/sum(IDX==12|IDX==8);
        k=k+1;
    end
end

eye1=[0,1,2,3];
k=1;
for i=1:length(eye1)
        ratio_dim_all(k)=sum(IDX'==8 & dyn_code(:,1)==eye1(i) )/sum(IDX==12|IDX==8);
        k=k+1;
end
groups={'0';'1';'2';'3';'00';'01';'02';'03';'10';'11';'12';'13';'20';'21';'22';'23';'30';'31';'32';'33'};
dimall=[ratio_dim_all,ratio_bindim_all];
[I,B]=sort(dimall,'descend');
figure;bar(1:1:20,dimall(B))
xticks(1:20)
set(gca,'xticklabel',groups(B))

%% for binocular dimmign responsive neurons

eye1=[0,1,2,3];
eye2=[0,1,2,3];
k=1;
for i=1:length(eye1)
    for j=1:length(eye2)
        ratio_bindim(k)=sum(IDX'==12 & dyn_code(:,1)==eye1(i) & dyn_code(:,2)==eye2(j))/sum(IDX==12);
        k=k+1;
    end
end
groups={'00';'01';'02';'03';'10';'11';'12';'13';'20';'21';'22';'23';'30';'31';'32';'33'};
[I,B]=sort(ratio_bindim,'descend');
figure;bar(1:1:16,ratio_bindim(B))
xticks(1:16)
set(gca,'xticklabel',groups(B))

%% for monocular dimmign responsive neurons

eye1=[0,1,2,3];
k=1;
for i=1:length(eye1)
        ratio_dim(k)=sum(IDX'==8 & dyn_code(:,1)==eye1(i) | IDX'==4 & dyn_code(:,3)==eye1(i))/sum(IDX==8|IDX==4);
        k=k+1;
end
groups={'0';'1';'2';'3'}
[I,B]=sort(ratio_dim,'descend');
figure;bar(1:1:4,ratio_dim(B))
xticks(1:4)
set(gca,'xticklabel',groups(B))

%% for looming
eye1=[0,1,2,3];
k=1;
for i=1:length(eye1)
        ratio_loom(k)=sum(IDX'==5 & dyn_code(:,2)==eye1(i) | IDX'==10 & dyn_code(:,1)==eye1(i))/sum(IDX==5|IDX==10);
        k=k+1;
end
groups={'0';'1';'2';'3'}
[I,B]=sort(ratio_loom,'descend');
figure;bar(1:1:4,ratio_loom(B))
xticks(1:4)
set(gca,'xticklabel',groups(B))
