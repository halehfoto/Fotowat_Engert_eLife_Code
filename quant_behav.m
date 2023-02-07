%This is a function that gets as input the flick data and then outputs
%response time, amplitude and direction.Created by Haleh Fotowat.
function [resp_time, Pmax, Pmin, PmaxminR, resp_amp, pm_amp, resp_dir,tail_dir_s,resp, pm_amp_isi,resp_dir_isi, resp_isi, resp_amp_isi, resp_time_isi,time_window] = quant_behav(flicks,time,stim_size_degrees,ISI)


%tail_dir_s=detrend((smooth(fillmissing(flicks,'spline'),10)),'constant');
tail_dir_s=fillmissing(flicks,'spline')-mean((fillmissing(flicks,'spline')));

%[pks_temp,locs_temp]=findpeaks(abs(tail_dir_s),'MinPeakHeight',max(max(abs(tail_dir_s))/2,20),'MinPeakProminence',10,'MinPeakWidth',5);
[pks_temp,locs_temp]=findpeaks(abs(tail_dir_s),'MinPeakHeight',15,'MinPeakProminence',10);

resp_time_temp=locs_temp/25;
time_window=time(min(find(stim_size_degrees==max(stim_size_degrees))))+5-time(find(stim_size_degrees>5,1));
index=find(resp_time_temp>time(find(stim_size_degrees>5,1)) & resp_time_temp<time(min(find(stim_size_degrees==max(stim_size_degrees))))+5);
dt=time(find(stim_size_degrees>5,1))-ISI;%time between the start of the stimulus and the starting point of accepting escapes
% index_isi=find(resp_time_temp<(ISI-dt) & resp_time_temp>(ISI-dt-time_window)) %look for tail flicks in the same size window prior to stim onset
index_isi=find(resp_time_temp<ISI & resp_time_temp>(ISI-time_window)); %look for tail flicks in the same size window prior to stim onset
locs=locs_temp(index);
pks=pks_temp(index);

locs_isi=locs_temp(index_isi);
pks_isi=pks_temp(index_isi);

if ~isempty(locs)
    resp=1;
    if tail_dir_s(locs(1))<0 
        resp_dir=1; %correct direction
    else
        resp_dir=-1;
    end
    
    resp_amp=tail_dir_s(locs(1));%amplitude of the first peak;
    pm_amp=tail_dir_s(locs(pks==max(pks)));%amplitude of the largest peak
    escape_thcross_1=find(abs(tail_dir_s(1:locs(1)))<10);%take the last point where it is less than 10 before escape starts
    if ~isempty(escape_thcross_1)
        resp_time=escape_thcross_1(end)/25;%timing of escape onset
    else
        resp_time=locs(1)/25;%response time this is used for clustering
        
    end
    
else
    resp=0;
    resp_dir=NaN;
    resp_time=NaN;
    resp_amp=0;
    pm_amp=0;
   
end
%do the same for ISI

if ~isempty(locs_isi)
    resp_isi=1;
    if tail_dir_s(locs_isi(1))<0 
        resp_dir_isi=1; %correct direction
    else
        resp_dir_isi=-1;
    end
    
    resp_amp_isi=tail_dir_s(locs_isi(1));%amplitude of the first peak;
    pm_amp_isi=tail_dir_s(locs_isi(pks_isi==max(pks_isi)));%amplitude of the largest peak
    escape_thcross_1_isi=find(abs(tail_dir_s(1:locs_isi(1)))<10);%take the last point where it is less than 10 before escape starts
    if ~isempty(escape_thcross_1_isi)
        resp_time_isi=escape_thcross_1_isi(end)/25;%timing of escape onset
    else
        resp_time_isi=locs_isi(1)/25;%response time this is used for clustering
        
    end
    
else
    resp_isi=0;
    resp_dir_isi=NaN;
    resp_time_isi=NaN;
    resp_amp_isi=0;
    pm_amp_isi=0;
   
end


Pmax=max(tail_dir_s(25*ISI+1:end));
Pmin=min(tail_dir_s(25*ISI+1:end));
PmaxminR=min(abs(Pmax),abs(Pmin))/max(abs(Pmax),abs(Pmin));
end