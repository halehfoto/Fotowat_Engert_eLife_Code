%this is a program that generates trial averages and forms clusters based
%on that/Choose folder : DMDLBRDLCB_allbrainpoolPL. Created by Haleh Fotowat.

close all
clearvars
path=uigetdir;
cd(path)

%read the pooled clustering data
load('DMDLBRDLCB_082628_091112_HA.mat')

RegData=csvread('DMDLBRDLCB_celltypes_registered.csv');
fish_n=RegData(:,1)
[NUM,TXT,RAW]=xlsread('DMDLBRDLCB_celltypes_registered_distribution.xlsx');
IDX=NUM(:,7);
ISI=40;
per2p=1.06295;
stim_type_cnt=5;
repse=5;
fps=69;
f0=round(ISI/(2*per2p));
nframes=size(F_dff_pool_all,2)
framerate=1/per2p;

time=1/framerate:1/framerate:nframes/framerate;

Length_trial_type=(size(F_dff_pool_all,2)/stim_type_cnt);
filenames=dir('*ha_allcells.mat');
%read the non-smooth F_dff
F_dff_pool_all_ns=[];
for i=1:length(filenames)
    load(filenames(i).name,'Keep_Units_All','F_dff','stim_size_degrees');
    F_dff_keep_All_ns{i}=(F_dff(Keep_Units_All{1},1:nframes));
    for ii=1:length(Keep_Units_All{1})
            F_dff_keep_ns_zs{i}(ii,:)=zscore(F_dff_keep_All_ns{i}(ii,:));
    end
    F_dff_pool_all_ns=[F_dff_pool_all_ns;F_dff_keep_ns_zs{i}];
end

[PeakAmps,PeakInd,PeakIndrelstim,expfit_peakt,confint_expt,expfit_signed,expfit_signed_bf] = CalcPeakResp(repse,stim_type_cnt,F_dff_pool_all,fps,ISI,per2p,time);
%for each cell calculate calculate the dimming and dark looming exponent

save('Peaks_EXPFits_DMDLBRDLCB_v2.mat','PeakAmps','PeakInd','PeakIndrelstim','expfit_peakt','confint_expt','expfit_signed','expfit_signed_bf','stim_size_degrees')

