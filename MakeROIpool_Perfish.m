%this is a program to take the plane number corresponding to the
%experimental plane and scale factor for all the fish in the pool and also the final pool
%cluster data and generate CSV files to be transformed into zbrain coordinates choose folder with relevant pool data
%Created by Haleh Fotowat
clearvars
close all
path=uigetdir;
cd(path)
filenames=dir('*AllBrains_procdata_ha_allcells.mat');
load('scaling information.mat')
for i=1:length(filenames)
    load(filenames(i).name,'A_keep','Keep_Units_All')
    units=Keep_Units_All{1};
    unitlocs=A_keep;
    for j=1:length(units)
        temp_cell=full(reshape(unitlocs(:,units(j)),1024,1034));
        temp_cell2=rot90(temp_cell,-1);
        props = regionprops(true(size(temp_cell2)), temp_cell2, 'WeightedCentroid');
        centroids{i}(j,:)=[props.WeightedCentroid(1),props.WeightedCentroid(2),iminfo(i,1)];
    end
    centroids{i}(:,1:2)=centroids{i}(:,1:2)*iminfo(i,2); %scale x and y by the scaling factor
    dlmwrite(strcat(filenames(i).name(1:12),'_Allbrain_ROIs_allcells.csv'),centroids{i},' ');
end

    
