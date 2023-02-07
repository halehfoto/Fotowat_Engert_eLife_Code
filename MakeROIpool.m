%this is a program to take the plane number corresponding to the
%experimental plane and scale factor for all the fish in the pool and also the final pool
%cluster data and generate CSV files to be transformed into zbrain coordinates choose folder DMDLBRDLCB_allbrainspoolPL
%created by Haleh Fotowat
clearvars
close all
path=uigetdir;
cd(path)

filenames=dir('*DFF.mat');
load('DMDLBRDLCB_082628_091112_HA.mat')
iminfo=[29	0.491
28	0.491
28	0.491
25	0.491
27	0.491
29	0.491
25	0.491
25	0.491
30	0.491
24	0.526
26	0.491
14  0.491  
16  0.491
16  0.491
14  0.491
];
for i=1:size(iminfo,1)
    temp1=find(fn'==i);
    units=cell_locs_all(temp1);
    unitlocs=Cells{i};
    for j=1:length(units)
        temp_cell=full(reshape(unitlocs(:,units(j)),1024,1034));
        temp_cell2=rot90(temp_cell,-1);
        props = regionprops(true(size(temp_cell2)), temp_cell2, 'WeightedCentroid');
        centroids{i}(j,:)=[props.WeightedCentroid(1),props.WeightedCentroid(2),iminfo(i,1)];
    end
    centroids{i}(:,1:2)=centroids{i}(:,1:2)*iminfo(i,2); %scale x and y by the scaling factor
    dlmwrite(strcat(filenames(i).name(1:12),'_Allbrain_ROIs_allcells.csv'),centroids{i},' ');
end

    
