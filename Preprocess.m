%this program reads the stimulus, behavior and 2p data and generates
%neuronal clusters for one plane recordings. Created by Haleh Fotowat
close all
clearvars
Reg=input('please enter 1 for analysis of hindbrain, 2 for midbrain and 3 for forebrain 0 for all brain and -1 for neuropil: ');
nCh=input('enter the number of 2p channels:');
taildata=input('Enter 1 if tail track data is available and zero if not:');
if nCh>1
    SO=input('enter the stim count after which the stimulus was turned off :');
end

if Reg==1
    Regname='HindBrains';
elseif Reg==2
    Regname='MidBrains';
elseif Reg==3
    Regname='ForeBrains';
elseif Reg==0
    Regname='AllBrains';
elseif Reg==-1
    Regname='NPBrains';
elseif Reg==4
    Regname='TL';
end
path=uigetdir;
temp=str2num(path(length(path)-1:length(path)));
if temp<10
    fishn=strcat('0',num2str(temp));
    b=0;
else
    fishn=num2str(temp);
    b=1;
end
filename=strcat(path(length(path)-(14+b):length(path)-(7+b)),'_f',fishn)

cd(path)
%% load stimulus file
% Indi=input('please enter the indicator used:');
% per2p=input('please enter imaging time per frame:');
Indi='gcamp7f';
per2p=1.06295;
switch Indi
    case 'gcamp7f'
        KernelWidth=0.962;
    case 'gcamp6s'
        KernelWidth=0.962;
    otherwise 
        KernelWidth=0.962;
end
cd stimulus
temp=dir('*stimulus_data.mat');
if temp(1).name(1)=='.'
    stimfname=temp(2);
else
    stimfname=temp(1);
end
    
[stim_size_degrees,stim_regress,isi_regress,time,nframes,fps,fpp,reps,repse,kernel,Stim_type,ISI,np]=Load_stim(per2p,stimfname,KernelWidth);
if exist('SO','var')
    SO_time=fps*SO*per2p;
    lightsoff_index_s=find(time>=SO_time,1);
else
    SO=NaN;
    lightsoff_index_s=length(stim_regress);
end
stim_regress(lightsoff_index_s:end)=median(stim_regress);
stim_size_degrees(lightsoff_index_s:end)=median(stim_size_degrees);
%generate a new time series that is one during stimulus + 10 seconds after
mask=zeros(size(stim_regress));
mask(stim_regress>median(stim_regress)+1)=1;
figure
plot(mask)
hold on
% plot(stim_regress)
for i=1:reps
    if (i*fps+ceil(5/per2p))<=fpp
        mask(i*fps:(i*fps+ceil(10/per2p)))=1;
    else
        mask(i*fps:fpp)=1;
    end
end
figure;plot(mask);hold on;plot(stim_regress/max(stim_regress)) 

if ~isnan(SO)
    stim_type_cnt=(SO-np)/repse;
else
    stim_type_cnt=(reps-np)/repse;
end
temp=[];
for i=1:stim_type_cnt
    if stim_type_cnt>1
        if cell2mat(Stim_type(:,np+(i-1)*repse+1))=='D'
            stmtp='DL'
            
        else
            stmtp=cell2mat(Stim_type(:,np+(i-1)*repse+1));
        end
        temp=strcat(temp,stmtp);
        
    else
        temp=Stim_type;
    end
end

if np>0
    temp=strcat('P',temp);
end

Exp_seq=temp(1:stim_type_cnt*2);
%% load behavior data
if taildata==1
    cd ../behavior
    %get behavior filenames
    temp=dir('*baseline.mat');
    if temp(1).name(1)=='.'
        behavfname=temp(2);
    else
        behavfname=temp(1);
    end
    
    flick_fill=[];
    flick_interp=[];
    vidrate=25;%video recording frame rate
    %calculate the extra number of frames that is acquired for each 2p plane
    %number of frames per stimuus
    tend2p=fps*per2p;
    peak_thresh=30;
    peak_dist=5;%200 ms
    peak_prom=10;
    for i=1:length(behavfname)
        load(behavfname(i).name)
        extra_frames_2p=round(5*per2p*25);
        ttail=1/vidrate:1/vidrate:(reps)*(length(flick(1,:))+extra_frames_2p)/vidrate;%time series for behavioral data
        if exist('SO','var')
            SO_time=fps*SO*per2p;
            lightsoff_index_b=find(ttail>=(SO_time-10),1);
        else
            SO=NaN;
            lightsoff_index_b=length(ttail);
        end
        flickpp_temp=[];
        for j=1:(reps)
            flickpp_temp=[flickpp_temp,[flick(i*j,:),flick(i*j,end)*ones(1,extra_frames_2p)]];
            %find out the response time and amplitude
            [amp,locp]=findpeaks(detrend((flick(i*j,:))),'MinPeakHeight',peak_thresh,'MinPeakDistance',peak_dist,'MinPeakProminence',peak_prom);
            [amp,locn]=findpeaks(detrend((-1*flick(i*j,:))),'MinPeakHeight',peak_thresh,'MinPeakDistance',peak_dist,'MinPeakProminence',peak_prom);
            tail_flick_idex{i,j}=[locp,locn];
            tail_flick_amp{i,j}=flick(i*j,tail_flick_idex{i,j});
        end
        flickpp_temp_fill=fillmissing(flickpp_temp,'spline');
        flickpp(i,:)=flickpp_temp_fill;
        flick_corrected(i,:)=flickpp(i,:);
        flick_corrected(i,lightsoff_index_b:end)=flick_corrected(i,lightsoff_index_b);
        flick_corrected(i,:)=flick_corrected(i,:)-mean(flick_corrected(i,:));
        
        temp=conv(flick_corrected,kernel,'full');
        motor_regress(i,:)=temp(1:length(flickpp_temp_fill));
        motor_regress2(i,:)=flickpp_temp_fill;
    end
    
    if length(behavfname)==1
        ztail=highpass(motor_regress(1,:),0.005);
        [pksb,locsb]=findpeaks(abs(ztail),'MinPeakProminence',30);
    end
    figure;plot(ttail,ztail);hold on;plot(ttail(locsb),ztail(locsb),'r*')
    hold on;plot(time,stim_regress)
    
    %% form a vector with stim type and behavior time. NaN if no response, first flick time if response
kkk=1;
resp_times=ttail(locsb);
resp_amps=ztail(locsb);
resp_time_series=zeros(size(ztail));
resp_time_series(locsb)=resp_amps;
%convovle response times with a gaussian of width 20 s (average width of
%the 20 responses)
kltime=-5:0.01:5;
% dbeta=betapdf(kltime,10,1)/max(betapdf(kltime,10,1));
dev=evpdf(kltime,0.2,0.5);
devg=gevpdf(kltime,0,1,0.5);

figure;plot(kltime,devg,'*')
% GW=gausswin(5);%1s wide gaussian window
temp=conv(resp_time_series,devg,'same');
behav_rate=temp(1:length(ttail));
figure(76);plot(ttail,behav_rate);hold on;plot(ttail,motor_regress(1,:));

%find response to the probe trial
    stim_time_trial=time(round(ISI/per2p):fps);
    resp_trial=find(resp_times>=stim_time_trial(1) &resp_times<=stim_time_trial(end));
    if ~isempty(resp_trial)
        behav_stim_pb=length(resp_trial);
        peak_time_rel_stimon_pb=resp_times(resp_trial)-stim_time_trial(1);
        peak_amp_pb=resp_amps(resp_trial);
    else
        behav_stim_pb=length(resp_trial);
        peak_time_rel_stimon_pb=resp_times(resp_trial)-stim_time_trial(1);
        peak_amp_pb=resp_amps(resp_trial);
        
    end
    for k=1:stim_type_cnt*repse
        
        %determine the response time
        stim_time_trial=time(round(ISI/per2p)+(k-1)*(fps)+np*(fps)+1:k*(fps)+np*fps);
        resp_trial=find(resp_times>=stim_time_trial(1) &resp_times<=stim_time_trial(end));
        
        
        if ~isempty(resp_trial)
            behav_stim_mat(k)=length(resp_trial);
            peak_time_rel_stimon{k}=resp_times(resp_trial)-stim_time_trial(1);
            peak_amp{k}=resp_amps(resp_trial);
        else
            peak_time_rel_stimon{k}=NaN;
            peak_amp{k}=NaN;
            behav_stim_mat(k)=0;
        end
        %overall response during ISI
        resp_isi=find(resp_times<stim_time_trial(1) & resp_times>stim_time_trial(1)-ISI);
        if ~isempty(resp_isi)
            behav_isi_mat(k)=length(resp_isi);
        else
            behav_isi_mat(k)=0;
        end
    end
    %calculate response rate during ISI
    spont_rate=behav_isi_mat./(ISI);
    %probability of spontaneous response during stimulus presentation
    alpha=spont_rate*(time(fps)-ISI);
    resp_amp_corr=(behav_stim_mat)-alpha;
    resp_amp_corr_pb=(behav_stim_pb)-alpha;
    figure
    plot(1:1:length(resp_amp_corr),resp_amp_corr,'*')
    figure
    plot(behav_stim_mat,'*')
    
    %save(strcat('/Volumes/LaCie/2P-Analyzed-backup/CaImAn-MATLAB-master/behavior_data/',num2str(repse),'reps/ISI',num2str(ISI),'/',filename,'_behav.mat'),'resp_amp_corr_pb','behav_stim_pb','peak_time_rel_stimon_pb','peak_amp_pb','peak_amp','peak_time_rel_stimon','resp_amp_corr','behav_stim_mat','behav_isi_mat','Stim_type','SO','np','reps','ISI');
end
%% get 2p filenames
cd ../looming/analyzed
filenames=dir('2*DFF.mat');
brainfilename=dir(strcat('2*',Regname,'ROIs.mat'));
% filenamev=dir('2*.h5');
% RawF=h5read(filenamev.name,'/mov');

%% for each plane now find units that either correlate with the stimulus or moto regressor 

xcorrthresh_s=0.3;
xcorrthresh_m=0.5;
xcorrthresh_t=0.7;

% xcorrthresh_s=0.1;
% xcorrthresh_m=0.5;
% xcorrthresh_t=0.1;

% resx=1024;
% resy=1034;
% first find traces that are correlated with stimulus or motor regressors
%take stimulus size directly, as the stimulus is so long in this case it
%doesnt make a difference if we use a regressor or the stimulus itself for
%correlations
stim_regress=stim_size_degrees;
for i=1:length(filenames)
    load(filenames(i).name,'F_dff','A_keep','ind_cnn','ind_corr','ind_exc','fitness','rval_space','value');
%     keep = (ind_corr | ind_cnn) & ind_exc;
%     fitness_keep=fitness(keep);
%     rval_space_keep=rval_space(keep);
%     value_keep=value(keep);
    
    load(brainfilename(i).name)
    [resx,resy]=size(brain);

    if Reg==1
        ROI=ROI1;
    elseif Reg==2
        ROI=ROI2;
    elseif Reg==3
        ROI=ROI3;
    elseif Reg==-1
        ROI=ROI0.Vertices;
    end
    F_dff_s=movmean(F_dff,10,2);%smooth F_dff
    [a,b]=size(A_keep);
    Keep_Units{i}=[];
    for j=1:b
        temp_cell=reshape(A_keep(:,j),resx,resy);
        [r,c]=find(temp_cell>0);
        if sum(inpolygon(c,r,ROI(:,1),ROI(:,2)))>0
            Keep_Units{i}=[Keep_Units{i},j];
        end
    end
    %find units that are correlated with any stimulus type

        
    index_temp=[];
    for j=1:stim_type_cnt
        [cs,lag]=corrcoef(([stim_regress((np-repse)*fps+j*repse*fps+1:np*fps+j*repse*fps);(F_dff(Keep_Units{i},(np-repse)*fps+j*repse*fps+1:np*fps+j*repse*fps))])');
        temp=find(abs(cs(1,:))>=xcorrthresh_s);
        index_temp=[index_temp,temp(2:end)-1];
        [ci,lag]=corrcoef(([isi_regress((np-repse)*fps+j*repse*fps+1:np*fps+j*repse*fps);(F_dff(Keep_Units{i},(np-repse)*fps+j*repse*fps+1:np*fps+j*repse*fps))])');
        temp=find(abs(ci(1,:))>=xcorrthresh_s);
        index_temp=[index_temp,temp(2:end)-1];

    end
    %stimulus responding neurons
    index_s=unique(index_temp);
    %units that show correlation with time
    [ct,lags]=corrcoef(([time(1:fpp).*stim_regress(1:fpp);((F_dff_s(Keep_Units{i},1:fpp)))])');
    temp=find(abs(ct(1,:))>=xcorrthresh_t);
    index_t=temp(2:end)-1;
    %index_temp=[index_temp,temp(2:end)-1]; %don't add these, time by
    %itself doesnt matter.
    if taildata==1
        %find motor responsive units
        temp_m=[];
        Motor=interp1(ttail,behav_rate,time,'spline');
        figure;plot(time,Motor);hold on;plot(ttail,motor_regress,'.')
        Ml=Motor;
        Ml(Motor>0)=0;
        plot(time,Ml,'b')
        Mr=Motor;
        Mr(Motor<=0)=0;
        plot(time,Mr,'r')
        index_m=[];
        index_m_lag=[];
        for j=1:size(F_dff_s(Keep_Units{i},:),1)
            %units that correlate with both tail movement direction
            [cmb,lagb]=xcorr(abs(Motor(1:lightsoff_index_s)),(F_dff_s(Keep_Units{i}(j),1:lightsoff_index_s)),20,'coeff');
            [cml,lagl]=xcorr(abs(Ml(1:lightsoff_index_s)),(F_dff_s(Keep_Units{i}(j),1:lightsoff_index_s)),20,'coeff');
            [cmr,lagr]=xcorr(abs(Mr(1:lightsoff_index_s)),(F_dff_s(Keep_Units{i}(j),1:lightsoff_index_s)),20,'coeff');
            if max(cmb)>=xcorrthresh_m | max(cml)>=xcorrthresh_m | max(cmr)>=xcorrthresh_m
                index_m=[index_m,j];
                [ab,bb]=max(cmb);
                [al,bl]=max(cml);
                [ar,br]=max(cmr);
                temp1=[ab,al,ar];
                temp2=[lagb(bb),lagl(bl),lagr(br)];
                
                [a,b]=max(temp1);
                index_m_lag=[index_m_lag,temp2(b)];%positive lag indicates that Motor peak comes after calcium peak
            end
        end
    else
        index_m=[];
    end
    
    index_temp=[index_temp,index_m];

    index_all=unique(index_temp);
    
    index_sm=find(ismember(index_m,index_s));
    index_mm=find(~ismember(index_m,index_s));
    Keep_Units_SM{i}=Keep_Units{i}(index_m(index_sm))
    Keep_Units_M{i}=Keep_Units{i}(index_m(index_mm))
    F_dff_keep_SM_nz{i}=(F_dff_s(Keep_Units_SM{i},1:nframes));
    F_dff_keep_M_nz{i}=(F_dff_s(Keep_Units_M{i},1:nframes));
    Keep_Units_S{i}=Keep_Units{i}(index_s)
    F_dff_keep_S_nz{i}=(F_dff_s(Keep_Units_S{i},1:nframes));
    Keep_Units_T{i}=Keep_Units{i}(index_t)
    F_dff_keep_T_nz{i}=(F_dff_s(Keep_Units_T{i},1:nframes));
    
    Keep_Units_All{i}=Keep_Units{i}(index_all);
    F_dff_keep_All_nz{i}=(F_dff_s(Keep_Units_All{i},1:nframes));
    
    for ii=1:length(index_s)
        F_dff_keep_S{i}(ii,:)=zscore(smooth(F_dff_keep_S_nz{i}(ii,:)));
    end
    for ii=1:(length(index_sm))
        F_dff_keep_SM{i}(ii,:)=zscore(smooth(F_dff_keep_SM_nz{i}(ii,:)));
    end
    for ii=1:length(index_mm)
            F_dff_keep_M{i}(ii,:)=zscore(smooth(F_dff_keep_M_nz{i}(ii,:)));
    end
    for ii=1:length(index_t)
            F_dff_keep_T{i}(ii,:)=zscore(smooth(F_dff_keep_T_nz{i}(ii,:)));
    end

    for ii=1:length(index_all)
            F_dff_keep{i}(ii,:)=zscore(smooth(F_dff_keep_All_nz{i}(ii,:)));
    end


    
end
close all
%%calculate the raw fluorescence
%% c-=luster using the whole waveformF_dff_pool=cell2mat(F_dff_keep');
% for i=1:length(Keep_Units_All{1})
%     temp_cell=reshape(A_keep(:,Keep_Units_All{1}(i)),resx,resy);
%     [pixx,pixy,MF]=find(temp_cell);
%     for k=1:size(F_dff_pool,2)
%         F_raw(i,k)=mean2(RawF(pixx,pixy,k));
%     end
% end
F_dff_pool=cell2mat(F_dff_keep');
F_dff_pool_nz=cell2mat(F_dff_keep_All_nz');
F_dff_pool_hp=F_dff_pool;
outfile=[];
for i=1:stim_type_cnt
    temp=Stim_type{(i-1)*repse+1+np};
    outfile=strcat(outfile,temp)
    if strmatch(temp,'TM')
        color_bar{i}='y';
    elseif strmatch(temp,'DM')
        color_bar{i}='k';
    elseif strmatch(temp,'DL')
        color_bar{i}='r';
    elseif strmatch(temp,'BR')
        color_bar{i}='c'
    elseif strmatch(temp,'CB')
        color_bar{i}='b';
    elseif strmatch(temp,'QF')
        color_bar{i}='m';
    elseif strmatch(temp,'DF')
        color_bar{i}='g';
    elseif strmatch(temp,'LR')
        color_bar{i}='r';
    elseif strmatch(temp,'LL')
        color_bar{i}='b';
    elseif strmatch(temp,'CR')
         color_bar{i}='m';
    elseif strmatch(temp,'CL')
         color_bar{i}='c';

    end
        
end
save(strcat(filename,'_',Regname,'_procdata_ha_allcells.mat'));
save(strcat('/Volumes/One Touch/2P-Analyzed-backup/CaImAn-MATLAB-master/PoolData/UniversalAnalysisHA/HC/DFDLDMCBBR/ISI40/',filename,'_',Regname,'_procdata_ha_allcells.mat'));

