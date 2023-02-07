%this is a program to plot the ROIs on the reference brain
%read the reference brain. Created by Haleh Fotowat.
close all
clearvars
path=uigetdir;
cd(path)
[X,meta]=nrrdread('Elavl3-H2BRFP.nrrd');
load('MaskDatabase.mat')


%read the pooled clustering data
load('DMDLBRDLCB_082628_091112_HA.mat')
fps=69;
ISI=40;
per2p=1.06295;
nframes=size(F_dff_pool_all,2)
framerate=1/per2p;
time=1/framerate:1/framerate:nframes/framerate;
RegData=csvread('DMDLBRDLCB_celltypes_registered_v2.csv');
fish_n=RegData(:,1)
[NUM,TXT,RAW]=xlsread('DMDLBRDLCB_celltypes_registered_v2_Hierarical_allnodes.xlsx');
colors=hsv(ccount);
%read the transformed ROIs

Y=std(single(X(:,:,:)),[],3);
Z=std(single(X(:,:,:)),[],1);
temp=(squeeze(Z))';
Zint=interp2(repmat(1:1:size(temp,1),size(temp,2),1),(repmat((1:1:size(temp,2))',1,size(temp,1))),temp',repmat(0.78/2:0.78/2:size(temp,1),size(temp,2),1),(repmat((1:1:size(temp,2))',1,length(0.78/2:0.78/2:size(temp,1)))));
Clusts=unique(RegData(:,5));
IDX=RegData(:,5);
EI=RegData(:,7);
%plot exciatatory vs inhibitory cells on the brain
ROI_pix=(NUM(:,1:3)./[0.798,0.798,0.798]);
%flag cells that are registered correctly and are associated to a region
flag=ones(size(IDX));
for i=1:(length(IDX)-1)
    if strcmp(TXT{i+1,8},'')
        flag(i)=0;
    end
end


for i=1:length(Clusts)
    figure(2*i)
    %figure(1)
    imagesc(Y);
    colormap gray
    hold on
    axis equal
    line([0 1400], [310.5,310.5])%approximate midline
    figure(2*i+1)
    %figure(2)
%     imagesc((squeeze(Z))');
    imagesc(Zint');

    colormap gray
    hold on
    set(gca,'YDir','normal')
    axis equal
    axis tight
end
Zplanes=unique(round(ROI_pix(:,3)));   
colors=jet(length(Clusts));
%%
%claculate the binocularity index
for i=1:length(Clusts)
    index=find(IDX==i-1 & flag==1 &fish_n~=10);
    %calculate the binocularity index
    BI(i)=(sum(ROI_pix(index,2)<310.5)-sum(ROI_pix(index,2)>310.5))/(sum(ROI_pix(index,2)<310.5)+sum(ROI_pix(index,2)>310.5));
end

%plot cells from each cluster
colormap1=1-jet
for i=1:length(Clusts)
    
    index=find(IDX==i & flag==1 &fish_n~=10);
    %calculate the binocularity index
    %for k=1:length(index)
        figure(2*i)
        plot(ROI_pix(index,1),ROI_pix(index,2),'LineStyle','none','Color',colormap1(17*Clusts(i)+1,:),'Marker','.','Markersize',6);
        title(strcat('Cluster number', num2str(i)))
        figure(2*i+1)
        plot(ROI_pix(index,1),ROI_pix(index,3),'LineStyle','none','Color',colormap1(17*Clusts(i)+1,:),'Marker','.','Markersize',6);
        title(strcat('Cluster number', num2str(i)))

%         figure(2)
%         plot(ROI_pix(index(k),1),ROI_pix(index(k),2),'Color','w','Marker','.','Markersize',6);
%         figure(3)
%         plot(ROI_pix(index(k),1),ROI_pix(index(k),3),'Color','w','Marker','.','Markersize',6);

    %end
%     index=find(IDX==i & flag==1 &fish_n~=10);
%     for k=1:length(index)
% %         figure(2*i)
% %         plot(ROI_pix(index(k),1),ROI_pix(index(k),2),'Color',colormap1(17*Clusts(i)+1,:),'Marker','.','Markersize',6);
% %         figure(2*i+1)
% %         plot(ROI_pix(index(k),1),ROI_pix(index(k),3),'Color',colormap1(17*Clusts(i)+1,:),'Marker','.','Markersize',6);
%     
%         figure(2*i)
%         plot(ROI_pix(index(k),1),ROI_pix(index(k),2),'Color','w','Marker','.','Markersize',6);
%         figure(2*i+1)
%         plot(ROI_pix(index(k),1),ROI_pix(index(k),3),'Color','w','Marker','.','Markersize',6);
% 
%     end
% 
    
    
    
end
figure(2)
imshow(Y',[]);hold on
temp_cell=reshape(full(MaskDatabaseOutlines(:,107)),[size(Y,2),size(Y,1),size(Z,3)]);
cell_loc_ind=find(temp_cell>0);
[loc_x,loc_y]=ind2sub(size(temp_cell),cell_loc_ind);
contour((max(temp_cell,[],3)),1,'LineColor','r','Linewidth',1);

%plot looming and dimming sensitive clusters together
% n=gcf;
% figure(n.Number+1)
% imagesc(Y);
% colormap gray
% hold on
% axis equal
% figure(n.Number+2)
% imagesc(Zint');
% colormap gray
% hold on
% set(gca,'YDir','normal')
% axis equal
% axis tight
figure(3)
imagesc(Y);
colormap gray
hold on
axis equal

plot(ROI_pix(flag==1 &fish_n~=10,1),ROI_pix(flag==1 &fish_n~=10,2),'LineStyle','none','Color','w','Marker','.','Markersize',4);

figure(4)
imagesc(Zint');
colormap gray
hold on
set(gca,'YDir','normal')
axis equal
axis tight
plot(ROI_pix(flag==1 &fish_n~=10,1),ROI_pix(flag==1 &fish_n~=10,3),'LineStyle','none','Color','w','Marker','.','Markersize',4);
figure(5)
imagesc(Y);
colormap gray
hold on
axis equal

plot(ROI_pix(flag==1 &fish_n~=10,1),ROI_pix(flag==1 &fish_n~=10,2),'LineStyle','none','Color','w','Marker','.','Markersize',4);

figure(6)
imagesc(Zint');
colormap gray
hold on
set(gca,'YDir','normal')
axis equal
axis tight
plot(ROI_pix(flag==1 &fish_n~=10,1),ROI_pix(flag==1 &fish_n~=10,3),'LineStyle','none','Color','w','Marker','.','Markersize',4);

for i=[5,12]
    
    index=find(IDX==i & flag==1 &fish_n~=10);

    if i==5
            figure(3)
            plot(ROI_pix(index,1),ROI_pix(index,2),'LineStyle','none','Color',colormap1(17*Clusts(i)+1,:),'Marker','.','Markersize',8);
    %         figure(n.Number+2)
            figure(4)
            plot(ROI_pix(index,1),ROI_pix(index,3),'LineStyle','none','Color',colormap1(17*Clusts(i)+1,:),'Marker','.','Markersize',8);
    elseif i==12
            figure(5)
            plot(ROI_pix(index,1),ROI_pix(index,2),'LineStyle','none','Color',colormap1(17*Clusts(i)+1,:),'Marker','.','Markersize',8);
    %         figure(n.Number+2)
            figure(6)
            plot(ROI_pix(index,1),ROI_pix(index,3),'LineStyle','none','Color',colormap1(17*Clusts(i)+1,:),'Marker','.','Markersize',8);

    end
    
end
%% Plot I cells
figure
imagesc(Y);
colormap gray
hold on
axis equal
plot(ROI_pix(flag==1 &fish_n~=10 &fish_n~=10,1),ROI_pix(flag==1 &fish_n~=10,2),'LineStyle','none','Color','w','Marker','.','Markersize',4);
plot(ROI_pix(flag==1 &fish_n~=10&EI==-1,1),ROI_pix(flag==1 &fish_n~=10&EI==-1,2),'LineStyle','none','Color','y','Marker','.','Markersize',8);

figure
imagesc(Zint');
colormap gray
hold on
set(gca,'YDir','normal')
axis equal
axis tight
plot(ROI_pix(flag==1 &fish_n~=10 &fish_n~=10,1),ROI_pix(flag==1 &fish_n~=10,3),'LineStyle','none','Color','w','Marker','.','Markersize',4);
plot(ROI_pix(flag==1 &fish_n~=10&EI==-1,1),ROI_pix(flag==1 &fish_n~=10&EI==-1,3),'LineStyle','none','Color','y','Marker','.','Markersize',8);

hold on
save('binocularity_index.mat','BI');