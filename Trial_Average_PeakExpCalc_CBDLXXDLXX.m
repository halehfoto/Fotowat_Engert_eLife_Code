%this is a program that generates trial averages and forms clusters based
%on that. Data in TimeAnalysis/XXCBDLXXDLXX. Created by Haleh Fotowat.

close all
clearvars
path=uigetdir;
cd(path)

stim_type_cnt=4;

% Length_trial_type=(size(F_dff_pool_all,2)/stim_type_cnt);
filenames=dir('*ha_allcells.mat');
%read the non-smooth F_dff
F_dff_pool_all=[];
for i=1:length(filenames)
    load(filenames(i).name,'repse','F_dff_pool','fps','ISI','per2p','time','stim_size_degrees');
    nframes=fps*stim_type_cnt*repse;
    stim_size_degrees=stim_size_degrees(1:nframes);
    F_dff_pool_all=[F_dff_pool_all;F_dff_pool(:,1:nframes)];
end

[PeakAmps,PeakInd,PeakIndrelstim,expfit_peakt,confint_expt,expfit_signed,expfit_signed_bf] = CalcPeakResp(repse,stim_type_cnt,F_dff_pool_all,fps,ISI,per2p,time);
%for each cell calculate calculate the dimming and dark looming exponent

save('Peaks_EXPFits_CBDLXXDLXX_v2.mat','PeakAmps','PeakInd','PeakIndrelstim','expfit_peakt','confint_expt','expfit_signed','expfit_signed_bf','stim_size_degrees','repse','fps','ISI','per2p','time')

