%This is a function that takes as input the tail direction information and
%plots a point for each trial in the three dimensional space with max+,
%max- and response duration and the color depicting response time with
%lighter earlier in the trial. Created by Haleh Fotowat
close all
% clearvars
rep=1;
s=what;
fb_e=[];%fin bear duration for escape trials
fb_ne=[];%fin beat duration for no escape trials
ibi_e=[];%interbeat intervals for escape trials
ibi_ne=[];%interbeat interval for no escape trials
for i=1:length(s.mat)
        load(s.mat{i});
        s.mat{i}
        figure;plot(Pl{3,2,1})
        if i==17
            break
        end
end

for i=1:length(s.mat)
    load(s.mat{i});
    if length(s.mat{i})==43
        isi((i-1)*10+1:i*10)=str2num(s.mat{i}(34:35));
    else
        isi((i-1)*10+1:i*10)=str2num(s.mat{i}(35:36));
    end
    for monitor=1:2
        if monitor==side
            k=1;
        else
            k=2;
        end
        for trial=1:10
            tail_temp{(i-1)*10+trial,k,rep}=tail_dir{trial,monitor,rep}-mean(tail_dir{trial,monitor,rep}(1:240));%subtract bas%line
            M((i-1)*10+trial,k)=monitor;
            %for each trial extract four variables min, max, response duration
            %and response time
            %             if ~isnan(escape_frame(trial,monitor,rep))
            %                 resp_var((i-1)*10+trial,k,rep,1)=ttc_vid_save{trial+(k-1)*10}(escape_frame(trial,monitor,rep))/1000;
            %             else
            %                 resp_var((i-1)*10+trial,k,rep,1)=ttc_vid_save{trial+(k-1)*10}(end)/1000;
            %             end
            mm=max(tail_temp{(i-1)*10+trial,k,rep});
            mmm=min(tail_temp{(i-1)*10+trial,k,rep});
            resp_var((i-1)*10+trial,k,rep,2)=max((tail_temp{(i-1)*10+trial,k,rep}));%max tail deflection
            resp_var((i-1)*10+trial,k,rep,3)=min(tail_temp{(i-1)*10+trial,k,rep});%min tail deflection
            resp_var((i-1)*10+trial,k,rep,4)=min(abs(mm),abs(mmm))/max(abs(mm),abs(mmm));%min max absolute ratio, if ratio is closer to 1 it is bidirectional
            %plot(tail_temp{(i-1)*10+trial,k,rep})
            tail_dir_s=smooth(tail_temp{(i-1)*10+trial,k,rep},5);

%             [pks,locs]=findpeaks(abs(tail_dir_s),'MinPeakHeight',max(max(abs(tail_dir_s))/2,20),'MinPeakProminence',10);
            [pks,locs]=findpeaks(abs(tail_dir_s),'MinPeakHeight',15,'MinPeakProminence',10);
            Npeaks((i-1)*10+trial,k,rep)=length(locs);
            if ~isempty(locs)
                if tail_dir_s(locs(1))<0 && monitor==1
                    resp_var((i-1)*10+trial,k,rep,11)=1; %correct direction
                elseif tail_dir_s(locs(1))>0 && monitor==2
                    resp_var((i-1)*10+trial,k,rep,11)=1; %correct direction
                else
                    resp_var((i-1)*10+trial,k,rep,11)=0;
                end

                resp_var((i-1)*10+trial,k,rep,8)=tail_dir_s(locs(1));%amplitude of the first peak;
                resp_var((i-1)*10+trial,k,rep,9)=tail_dir_s(locs(pks==max(pks)));%amplitude of the largest peak
                
                escape_thcross_1=find(abs(tail_dir_s(1:locs(1)))<10);%take the last point where it is less than 10 before escape starts
                escape_thcross_2=find(abs(tail_dir_s(locs(end):end))>10);%take the last point where it is more than 10 after the last peak
                if ~isempty(escape_thcross_1)
                    resp_var((i-1)*10+trial,k,rep,5)=(locs(end)-escape_thcross_1(end))/120;%time between escape onset and the last peak
                    resp_var((i-1)*10+trial,k,rep,1)=escape_thcross_1(end)/120;%timing of escape onset
                    if (round(escape_thcross_1(end)/6)+1)<=length(ang_size_save{trial+(monitor-1)*10})
                        resp_var((i-1)*10+trial,k,rep,7)=ang_size_save{trial+(monitor-1)*10}((round(escape_thcross_1(end)/6)+1));
                    else
                        resp_var((i-1)*10+trial,k,rep,7)=ang_size_save{trial+(monitor-1)*10}(end);%if response happened after expansion
                    end
                else
                    resp_var((i-1)*10+trial,k,rep,1)=locs(1)/120;%response time this is used for clustering
                    resp_var((i-1)*10+trial,k,rep,5)=(locs(end)-locs(1))/120;
                    if (round(locs(1)/6)+1)<=length(ang_size_save{trial+(monitor-1)*10})
                        resp_var((i-1)*10+trial,k,rep,7)=ang_size_save{trial+(monitor-1)*10}(round(locs(1)/6)+1);
                    else
                        resp_var((i-1)*10+trial,k,rep,7)=ang_size_save{trial+(monitor-1)*10}(end);%if response happened after expansion
                    end
                    
                end
                resp_var((i-1)*10+trial,k,rep,10)=resp_var((i-1)*10+trial,k,rep,1)+(ttc_vid_save{trial+(monitor-1)*10}(1)/1000);%response time relative to collision

            else
                resp_var((i-1)*10+trial,k,rep,5)=0;%respose duration
                resp_var((i-1)*10+trial,k,rep,1)=-1;%response time =-1 if no response detected i.e. no peak this is used for clustering
                resp_var((i-1)*10+trial,k,rep,7)=NaN;
                resp_var((i-1)*10+trial,k,rep,8)=0;%amplitude of the first peak;
                resp_var((i-1)*10+trial,k,rep,9)=mm;%largest amplitude
                resp_var((i-1)*10+trial,k,rep,10)=NaN;%response time relative to collision 
                resp_var((i-1)*10+trial,k,rep,11)=NaN;
            end
            
            resp_var((i-1)*10+trial,k,rep,6)=trial;%trial number
            %calculate the fin-flick rate during the time the stimulus is
            %less than 10 degrees and the time in between 10 degrees and
            %maximum size
            %find the timing of the last fin movement, in escape trials,
            %make sure that this is taken as the last one before escape
            if resp_var((i-1)*10+trial,k,rep,1)~=-1 && ~isempty(escape_thcross_1)
                FinL=max(locs_fin_only{trial,monitor,rep}(locs_fin_only{trial,monitor,rep}<escape_thcross_1(end)));
            elseif resp_var((i-1)*10+trial,k,rep,1)~=-1 && isempty(escape_thcross_1)
                FinL=max(locs_fin_only{trial,monitor,rep}(locs_fin_only{trial,monitor,rep}<locs(1)));
            else
                FinL=max(locs_fin_only{trial,monitor,rep}(locs_fin_only{trial,monitor,rep}<find(ttc_vid_save{trial+(monitor-1)*10}>0,1)));
            end
%             if FinL
%             FinL
%             figure;plot(tail_dir_s);hold on;plot(Pl{trial,monitor,rep});plot(Pr{trial,monitor,rep});
%             
%             end
            if ~isempty(FinL)
                    resp_var((i-1)*10+trial,k,rep,13)=(FinL/120)+ttc_vid_save{trial+(monitor-1)*10}(1)/1000;%timing of the last fin movement relative to collision
                else
                    resp_var((i-1)*10+trial,k,rep,13)=NaN;
            end
            FinR1=sum(locs_fin_only{trial,monitor,rep}<find(ttc_vid_save{trial+(monitor-1)*10}>-2000,1))/(-2-ttc_vid_save{trial+(monitor-1)*10}(1)/1000);%rate per second
            FinR2=sum(locs_fin_only{trial,monitor,rep}>find(ttc_vid_save{trial+(monitor-1)*10}>-2000,1)& locs_fin_only{trial,monitor,rep}<find(ttc_vid_save{trial+(monitor-1)*10}>0,1))/2;%rate per second
            if FinR1>0 
                resp_var((i-1)*10+trial,k,rep,12)=FinR2/FinR1;%the ratio of fins before and after 2 s before collision
            elseif FinR1==0 && FinR2>0
                resp_var((i-1)*10+trial,k,rep,12)=10;%
            else
                resp_var((i-1)*10+trial,k,rep,12)=-10;
            end
            resp_var((i-1)*10+trial,k,rep,12)
            
%             cla
            %properties of finbeats
            if resp_var((i-1)*10+trial,k,rep,1)~=-1
                esctime=resp_var((i-1)*10+trial,k,rep,1)*120-120;
            else
                esctime=length(ttc_vid_save{trial});
            end
            temp_finr=Pr{trial,monitor,rep}(1:esctime);
            temp_finr=temp_finr-mean(temp_finr);
            temp_finl=Pl{trial,monitor,rep}(1:esctime);
            temp_finl=temp_finl-mean(temp_finl);
            if ~isempty(temp_finr) && ~isempty(temp_finl)
%                 pulsewidth(smooth(temp_finr),ttc_vid_save{trial}(1:esctime),'Polarity','Positive','StateLevels',[-1,4],'Tolerance',30);
                BD_2r=pulsewidth(smooth(temp_finr),ttc_vid_save{trial}(1:esctime),'Polarity','Positive','StateLevels',[-1,4],'Tolerance',30);
%                 pulsewidth(smooth(temp_finl),ttc_vid_save{trial}(1:esctime),'Polarity','Positive','StateLevels',[-1,4],'Tolerance',30);
                BD_2l=pulsewidth(smooth(temp_finl),ttc_vid_save{trial}(1:esctime),'Polarity','Positive','StateLevels',[-1,4],'Tolerance',30);
                pk_thresh=10;
                pk_prom=5;
                [pks_finr,locs_finr]=findpeaks(temp_finr,'MinPeakHeight',pk_thresh,'MinPeakProminence',pk_prom);%find all the peaks
                [pks_finl,locs_finl]=findpeaks(temp_finl,'MinPeakHeight',pk_thresh,'MinPeakProminence',pk_prom);%find all the peaks
                IbeatI{(i-1)*10+trial,k,rep}=1000*[diff(locs_finr),diff(locs_finl)]/120;%inter fin beat interval in ms
                beat_dur_dist{(i-1)*10+trial,k,rep}=[BD_2r;BD_2l];%fin beating duration (ms)
            else
                IbeatI{(i-1)*10+trial,k,rep}=[];
                beat_dur_dist{(i-1)*10+trial,k,rep}=[];
            end
            if resp_var((i-1)*10+trial,k,rep,1)~=-1 &k==1
                fb_e=[fb_e;beat_dur_dist{(i-1)*10+trial,k,rep}];%fin bear duration for escape trials
                ibi_e=[ibi_e,IbeatI{(i-1)*10+trial,k,rep}];%interbeat intervals for escape trials
            elseif resp_var((i-1)*10+trial,k,rep,1)==-1 &k==1

                fb_ne=[fb_ne;beat_dur_dist{(i-1)*10+trial,k,rep}];;%fin beat duration for no escape trials
                ibi_ne=[ibi_ne,IbeatI{(i-1)*10+trial,k,rep}];%interbeat interval for no escape trials
            end
        tail_trace{(i-1)*10+trial,k,rep}=tail_dir_s;
        if  ang_size_save{trial}(find(ttc_stim_save{trial}/1000>resp_var((i-1)*10+trial,k,rep,10),1))>10
            flag_es((i-1)*10+trial,k,rep)=1;%escape vs struggle
        else
             flag_es((i-1)*10+trial,k,rep)=-1;
        end
        end
%         figure(1)
%         plot(ttc_vid_save{trial+(monitor-1)*10}(1:length(tail_dir_s))/1000,abs(tail_dir_s))
%         hold on
%         plot(ttc_stim_save{trial+(monitor-1)*10}(1:length(ang_size_save{trial+(monitor-1)*10}))/1000,ang_size_save{trial+(monitor-1)*10})
%         resp_var((i-1)*10+trial,k,rep,7)
%         close all
    end
    i
end
%compare distribution of fin beat duration for escape and no escape trials
fb=[fb_ne;fb_e];
fbg=[ones(size(fb_ne));2*ones(size(fb_e))];
kruskalwallis(fb,fbg)
dfittool(ibi_e)
dfittool(ibi_ne)
%show the distribution of fin and tail times only for the first trial in
%each fish where I had both
for i=1:length(s.mat)
    
    index_1(i)=(i-1)*10+1;
end
index=find(~isnan(resp_var(index_1,1,1,13))& ~isnan(resp_var(index_1,1,1,10)));

dfittool(resp_var(index_1(index),1,1,13))
dfittool(resp_var(index_1(index),1,1,10))

dfittool(resp_var(~ismember(index_1,index),1,1,10))


%plot the timing of fin stopping for freeze trials vs escape trials
index_freeze_1=find(~isnan(resp_var(:,1,1,13))&isnan(resp_var(:,1,1,10)));
index_freeze_2=find(~isnan(resp_var(:,2,1,13))&isnan(resp_var(:,2,1,10)));

index_escape_1=find(~isnan(resp_var(:,1,1,13))&~isnan(resp_var(:,1,1,10)));
index_escape_2=find(~isnan(resp_var(:,2,1,13))&~isnan(resp_var(:,2,1,10)));

dfittool([resp_var(index_freeze_1,1,1,13)])
dfittool([resp_var(index_escape_1,1,1,13)])
dfittool([resp_var(index_escape_1,1,1,10)])
%compare the timing of fin tuck in escape and freeze trials, it happened
%significantly earlier in escape trials.

kruskalwallis([resp_var(index_freeze_1,1,1,13);resp_var(index_escape_1,1,1,13)],[ones(size(resp_var(index_freeze_1,1,1,13)));2*ones(size(resp_var(index_escape_1,1,1,13)))])


foldername = 'tap data';
         % folder where all the files are located.
filetype = 'mat'; % type of files to be processed
        % Types currently supported .tif/.tiff, .h5/.hdf5, .raw, .avi, and .mat files
files_temp = dir(fullfile(foldername,['*.',filetype]));   % list of filenames (will search all subdirectories)
for i=1:length(files_temp)
    load(fullfile(foldername,files_temp(i).name))
    if length(files_temp(i).name)==47
        isi_tap(i)=str2num(files_temp(i).name(34:35));%the isi belongs to the isi of the looming trials that came before
    else
        isi_tap(i)=str2num(files_temp(i).name(35:36));
    end
    
    tail_tempt{i}=tail_dir{1}-mean(tail_dir{1}(250:end));%subtract bas%line, which is the end of recording in the case of tap
%     mm=max(tail_tempt{i});
%     mmm=min(tail_tempt{i});
    tail_dir_s=tail_tempt{i};
    %tap stimulus comes on the 20th frame
    [pks,locs]=findpeaks(abs(tail_dir_s(20:end)),'MinPeakHeight',15,'MinPeakProminence',10);
    mm=max(tail_dir_s(20:end));
    mmm=min(tail_dir_s(20:end));
    if ~isempty(locs)
        resp_var_tap(i)=tail_dir_s((20+locs(1)));%amplitude of the first peak
        %tap was presented on the 20th frame
        index=find(abs(tail_dir_s(1:locs(1)))<10);%take the last point where it is less than 10 before escape starts
        if ~isempty(index) 
            resp_var_tapt(i)=1000*(index(end))/120;%time rel to tap in ms
        else
            resp_var_tapt(i)=1000*(locs(1))/120;
        end
    else
        resp_var_tap(i)=mean(tail_dir_s(250:end));
        resp_var_tapt(i)=NaN;
    end
    i
end

%plot the response amplitude vs trial number for the 1st and the second
%side
%find  largest peaks and use its absolute value

resp_amps_s1=abs(resp_var(:,1,1,9));
resp_amps_s2=abs(resp_var(:,2,1,9));

resp_amp_tap=abs(resp_var_tap)

resp_times_s1=(resp_var(:,1,1,10));
resp_times_s2=(resp_var(:,2,1,10));

resp_FinR_s1=(resp_var(:,1,1,12));
resp_FinR_s2=(resp_var(:,2,1,12));

resp_FinR_E_s1=(resp_var(resp_var(:,1,1,1)>-1,1,1,12));%trials with escape
resp_FinR_E_s2=(resp_var(resp_var(:,1,1,1)>-1,2,1,12));

resp_FinR_NE_s1=(resp_var(resp_var(:,1,1,1)==-1,1,1,12));%trials with escape
resp_FinR_NE_s2=(resp_var(resp_var(:,1,1,1)==-1,2,1,12));

isi(isi==3)=180;
ISI=[5,10,180];
NFISH=[sum(isi==5),sum(isi==10),sum(isi==180)]/10;
for i=1:length(ISI)
    for trials=1:10
        Ave_amps_s1(i,trials,1)=nanmean(resp_amps_s1(resp_var(:,1,1,6)==trials & (isi==ISI(i))'));
        Ave_amps_s1(i,trials,2)=nanstd(resp_amps_s1(resp_var(:,1,1,6)==trials & (isi==ISI(i))'));
        
        Ave_times_s1(i,trials,1)=nanmean(resp_times_s1(resp_var(:,1,1,6)==trials & (isi==ISI(i))'));
        Ave_times_s1(i,trials,2)=nanstd(resp_times_s1(resp_var(:,1,1,6)==trials & (isi==ISI(i))'));

        Ave_amps_s2(i,trials,1)=nanmean(resp_amps_s2(resp_var(:,2,1,6)==trials & (isi==ISI(i))'));
        Ave_amps_s2(i,trials,2)=nanstd(resp_amps_s2(resp_var(:,2,1,6)==trials & (isi==ISI(i))'));
        
        Ave_times_s2(i,trials,1)=nanmean(resp_times_s2(resp_var(:,2,1,6)==trials & (isi==ISI(i))'));
        Ave_times_s2(i,trials,2)=nanstd(resp_times_s2(resp_var(:,2,1,6)==trials & (isi==ISI(i))'));
       
        Ave_FinR_s1(i,trials,1)=nanmean(resp_FinR_s1(resp_var(:,1,1,6)==trials & (isi==ISI(i))'));
        Ave_FinR_s1(i,trials,2)=nanstd(resp_FinR_s1(resp_var(:,1,1,6)==trials & (isi==ISI(i))'));
      
        Ave_FinR_s2(i,trials,1)=nanmean(resp_FinR_s2(resp_var(:,2,1,6)==trials & (isi==ISI(i))'));
        Ave_FinR_s2(i,trials,2)=nanstd(resp_FinR_s2(resp_var(:,2,1,6)==trials & (isi==ISI(i))'));
        
        
        Ave_FinR_s1_allisi(trials,1)=nanmean(resp_FinR_s1(resp_var(:,1,1,6)==trials));
        Ave_FinR_s1_allisi(trials,2)=nanstd(resp_FinR_s1(resp_var(:,1,1,6)==trials));
      
        Ave_FinR_s2_allisi(trials,1)=nanmean(resp_FinR_s2(resp_var(:,2,1,6)==trials));
        Ave_FinR_s2_allisi(trials,2)=nanstd(resp_FinR_s2(resp_var(:,2,1,6)==trials));

        %Average trials with escape
        Ave_FinR_s1_allisi_E(trials,1)=nanmean(resp_FinR_s1(resp_var(:,1,1,6)==trials & resp_var(:,1,1,1)~=-1 ));
        Ave_FinR_s1_allisi_E(trials,2)=nanstd(resp_FinR_s1(resp_var(:,1,1,6)==trials & resp_var(:,1,1,1)~=-1 ));
      
        Ave_FinR_s2_allisi_E(trials,1)=nanmean(resp_FinR_s2(resp_var(:,2,1,6)==trials & resp_var(:,2,1,1)~=-1 ));
        Ave_FinR_s2_allisi_E(trials,2)=nanstd(resp_FinR_s2(resp_var(:,2,1,6)==trials & resp_var(:,2,1,1)~=-1));
        
        %Average trials with no escape
        Ave_FinR_s1_allisi_NE(trials,1)=nanmean(resp_FinR_s1(resp_var(:,1,1,6)==trials & resp_var(:,1,1,1)==-1));
        Ave_FinR_s1_allisi_NE(trials,2)=nanstd(resp_FinR_s1(resp_var(:,1,1,6)==trials & resp_var(:,1,1,1)==-1));
      
        Ave_FinR_s2_allisi_NE(trials,1)=nanmean(resp_FinR_s2(resp_var(:,2,1,6)==trials & resp_var(:,2,1,1)==-1 ));
        Ave_FinR_s2_allisi_NE(trials,2)=nanstd(resp_FinR_s2(resp_var(:,2,1,6)==trials & resp_var(:,2,1,1)==-1 ));
        
        
        pcnt_crct(i,trials,1)=nansum(resp_var(resp_var(:,1,1,6)==trials & (isi==ISI(i))',1,1,11)==1)/sum(~isnan(resp_var(resp_var(:,1,1,6)==trials & (isi==ISI(i))',1,1,11)));
        pcnt_crct(i,trials,2)=nansum(resp_var(resp_var(:,2,1,6)==trials & (isi==ISI(i))',2,1,11)==1)/sum(~isnan(resp_var(resp_var(:,2,1,6)==trials & (isi==ISI(i))',2,1,11)));
    end

    figure(1)
    
    subplot(3,1,i)
    title(strcat('isi=',num2str(ISI(i)),'s', '  nfish=', num2str(NFISH(i))))
    hold on
    plot(1:1:10,Ave_amps_s1(i,:,1),'r*');
    f=fit((1:1:10)',Ave_amps_s1(i,:,1)','exp1');
    cint=confint(f);
    fitbias1(i,:)=[f.a,cint(:,1)']
    fitexpon1(i,:)=[f.b,cint(:,2)']
    hl1=plot(f,(1:1:10)',Ave_amps_s1(i,:,1)')
    plot(1:1:10,Ave_amps_s1(i,:,1)+Ave_amps_s1(i,:,2)./sqrt(sum(isi==ISI(i))/10),'r:');
    plot(1:1:10,Ave_amps_s1(i,:,1)-Ave_amps_s1(i,:,2)./sqrt(sum(isi==ISI(i))/10),'r:');
    set([hl1(1) hl1(2)],'color','r')
    plot(11:1:20,Ave_amps_s2(i,:,1),'k*');
    f=fit((11:1:20)',Ave_amps_s2(i,:,1)','exp1');
    cint=confint(f);
    fitbias2(i,:)=[f.a,cint(:,1)']
    fitexpon2(i,:)=[f.b,cint(:,2)']

    hl2=plot(f,(11:1:20)',Ave_amps_s2(i,:,1)')
    set([hl2(1) hl2(2)],'color','k')

    plot(11:1:20,Ave_amps_s2(i,:,1)+Ave_amps_s2(i,:,2)./sqrt(sum(isi==ISI(i))/10),'k:');
    plot(11:1:20,Ave_amps_s2(i,:,1)-Ave_amps_s2(i,:,2)./sqrt(sum(isi==ISI(i))/10),'k:');
    legend off
    box off
    %line([0.95 10.05], [20,20])
    %plot the tap response
    %errorbar(21,nanmean(resp_amp_tap(isi_tap==ISI(i))),nanstd(resp_amp_tap(isi_tap==ISI(i)))/(sqrt(length(isi_tap==ISI(i)))))
    
    xlim([0.95 20.05])
    if i==3
        xlabel('trial number')
    end
    if i==2
        ylabel({'Average tail';'deviation from mid%line';'(pixels)'})
    end

    text(7,50,strcat('exp= ', num2str(fitexpon1(i,1),'%.2f')),'Color','r')
    text(17,50,strcat('exp= ', num2str(fitexpon2(i,1),'%.2f')),'Color','k')
    figure(2)
    subplot(3,1,i)
    title(strcat('isi=',num2str(ISI(i)),'s', '  nfish=', num2str(NFISH(i))))

    hold on
    plot(1:1:10,Ave_times_s1(i,:,1),'r');
    plot(1:1:10,Ave_times_s1(i,:,1)+Ave_times_s1(i,:,2)./sqrt(sum(isi==ISI(i))/10),'r:');
    plot(1:1:10,Ave_times_s1(i,:,1)-Ave_times_s1(i,:,2)./sqrt(sum(isi==ISI(i))/10),'r:');
    plot(11:1:20,Ave_times_s2(i,:,1),'k');
    plot(11:1:20,Ave_times_s2(i,:,1)+Ave_times_s2(i,:,2)./sqrt(sum(isi==ISI(i))/10),'k:');
    plot(11:1:20,Ave_times_s2(i,:,1)-Ave_times_s2(i,:,2)./sqrt(sum(isi==ISI(i))/10),'k:');
    ylim([round(ttc_vid_save{1}(1))/1000,0.5])
    line([0.95 20.05],[0,0])
    xlim([0.95 20.05])
    if i==3
        xlabel('trial number')
    end
    if i==2
    ylabel({'Response time';'rel to collision';'(s)'})
    end
%     subplot(4,1,3)
%     hold on
%     plot(1:1:10,Ave_times_s1(i,:,2),'r');
%     plot(1:1:10,Ave_times_s2(i,:,2),'k');
%     
%     p_SD(i)=kruskalwallis([Ave_times_s1(i,5:10,2),Ave_times_s2(i,5:10,2)],[ones(size(Ave_times_s1(i,5:10,2))),2*ones(size(Ave_times_s2(i,5:10,2)))],'off');
    figure(3)
    
    subplot(3,1,i)
    title(strcat('isi=',num2str(ISI(i)),'s', '  nfish=', num2str(NFISH(i))))
    hold on
    plot(1:1:10,pcnt_crct(i,:,1),'r*-')
    hold on
    plot(11:1:20,pcnt_crct(i,:,2),'k*-')
    ylim([0,1])
    line([0.95 20.05],[0.5,0.5])
    legend off
    box off
    if i==3
        xlabel('trial number')
    end
    if i==2
        ylabel({'correct dir';'probability'})
    end
    xlim([0.95 20.05])
        figure(4)
    subplot(3,1,i)
    hold on
    plot(1:1:10,Ave_FinR_s1(i,:,1),'g');
    plot(1:1:10,Ave_FinR_s1(i,:,1)+Ave_FinR_s1(i,:,2)./sqrt(sum(isi==ISI(i))/10),'g:');
    plot(1:1:10,Ave_FinR_s1(i,:,1)-Ave_FinR_s1(i,:,2)./sqrt(sum(isi==ISI(i))/10),'g:');
    
    plot(11:1:20,Ave_FinR_s2(i,:,1),'b');
    plot(11:1:20,Ave_FinR_s2(i,:,1)+Ave_FinR_s2(i,:,2)./sqrt(sum(isi==ISI(i))/10),'b:');
    plot(11:1:20,Ave_FinR_s2(i,:,1)-Ave_FinR_s2(i,:,2)./sqrt(sum(isi==ISI(i))/10),'b:');
    line([0.95 20.05],[0,0])


end
figure(1)
saveas(gcf,'resp_amp_perisi.pdf')
figure(2)
saveas(gcf,'resp_time_perisi.pdf')
figure(3)
saveas(gcf,'resp_corr_perisi.pdf')
figure(4)
saveas(gcf,'resp_fin_perisi.pdf')

%pool ISIs and pool the first three (F3) and last three responses  (L3)
F3S1=resp_amps_s1(resp_var(:,1,1,6)==1 | resp_var(:,1,1,6)==2 | resp_var(:,1,1,6)==3);
F3S2=resp_amps_s2(resp_var(:,2,1,6)==1 | resp_var(:,2,1,6)==2 | resp_var(:,2,1,6)==3);

L3S1=resp_amps_s1(resp_var(:,1,1,6)==8 | resp_var(:,1,1,6)==9 | resp_var(:,1,1,6)==10);
L3S2=resp_amps_s2(resp_var(:,2,1,6)==8 | resp_var(:,2,1,6)==9 | resp_var(:,2,1,6)==10);

TAP=resp_amp_tap;

COMP_AMP=[F3S1',L3S1',F3S2',L3S2',TAP];
Groups=[ones(size(F3S1')),2*ones(size(L3S1')),3*ones(size(F3S2')),4*ones(size(L3S2')),5*ones(size(TAP))];
[p,a,stat]=kruskalwallis(COMP_AMP,Groups);
figure
p=multcompare(stat)


figure
plot(1:1:10, Ave_FinR_s1_allisi(:,1),'r');
hold on
plot(1:1:10, Ave_FinR_s2_allisi(:,1),'k');


plot(1:1:10, Ave_FinR_s1_allisi_E(:,1),'b');

figure
plot(1:1:10, Ave_FinR_s2_allisi_NE(:,1),'k');
hold on
plot(1:1:10, Ave_FinR_s2_allisi_E(:,1),'r');

figure(1)
saveas(gcf,strcat('resp_amp',num2str(ISI(i)),'.pdf'))
saveas(gcf,strcat('resp_amp',num2str(ISI(i)),'.png'))
figure(2)
saveas(gcf,strcat('resp_time',num2str(ISI(i)),'.pdf'))
saveas(gcf,strcat('resp_time',num2str(ISI(i)),'.png'))
figure(3)
saveas(gcf,strcat('resp_dir',num2str(ISI(i)),'.pdf'))
saveas(gcf,strcat('resp_dir',num2str(ISI(i)),'.png'))

%find the 10 and 90 percent quantile of escape time data
X=[resp_var(:,1,1,10);resp_var(:,2,1,10)];
Q=quantile(X,[0.1,0.5,0.90]);

dfittool(X)
dfittool(X(X>Q(1)&X<Q(3)))

Ang=[resp_var(:,1,1,7);resp_var(:,2,1,7)]
dfittool(Ang(X>Q(1)&X<Q(3)));

M_reshape=reshape(M,[length(s.mat)*20,1]);
ngroups=5;
y=[resp_var(:,1,1,1:4);resp_var(:,2,1,1:4)];
y_trials=[resp_var(:,1,1,6);resp_var(:,2,1,6)];
idx=kmeans(squeeze(y),ngroups);
idx_m(1,:)=idx(1:10*length(s.mat));
y_trials_m(1,:)=y_trials(1:10*length(s.mat));
idx_m(2,:)=idx(1+10*length(s.mat):20*length(s.mat));
y_trials_m(2,:)=y_trials(1+10*length(s.mat):20*length(s.mat));
%Do a different gouping based on stimulation side left eye or right eye
y_trials_s(1,:)=y_trials(M_reshape==1);
y_trials_s(2,:)=y_trials(M_reshape==2);
idx_s(1,:)=idx(M_reshape==1);
idx_s(2,:)=idx(M_reshape==2);



%%
colors=jet(20);
%side 1 and side 2
CNT=0;
for i=1:length(s.mat)
    load(s.mat{i});
    for monitor=1:2
        if monitor==side
            k=1;
        else
            k=2;
        end
        
        for trial=1:10
            figure(40+k)
            %colors=colormap(hot(ngroups+3));
            %if abs(resp_var((i-1)*10+trial,k,rep,12))~=10
                plot3(resp_var((i-1)*10+trial,k,rep,2),resp_var((i-1)*10+trial,k,rep,3),resp_var((i-1)*10+trial,k,rep,4),'Color',colors(idx_m(k,(i-1)*10+trial)*4,:),'Marker','.','MarkerSize',10);
                hold on

                text(resp_var((i-1)*10+trial,k,rep,2),resp_var((i-1)*10+trial,k,rep,3),resp_var((i-1)*10+trial,k,rep,4),num2str(idx_m(k,(i-1)*10+trial)))
                grid on
                xlabel('Right flick')
                ylabel('Left flick')
                zlabel('bidirectionality index')
%                 h=colorbar;
%                 h.Ticks=1/(2*(ngroups+3)):1/(ngroups+3):(ngroups+2)/(ngroups+3)+1/(2*(ngroups+3))
%                 h.TickLabels={'1','2','3','4','5','6','7','8'};
                %ylim([0,2])
            %end
            
        end
    end
end
figure(41)
saveas(gcf,'ISIAll clusters side1.pdf')
figure(42)
saveas(gcf,'ISIAll clusters side2.pdf')
Group_names={'Largebi','noflick', 'left','right','swim'};
%do the same thing separating monitor one and two
[x1,y1]=find(M==1);
[x2,y2]=find(M==2);
figure(71)
colors=colormap(hot(ngroups+3));

for j=1:length(x1)
    plot3(resp_var(x1(j),y1(j),rep,2),resp_var(x1(j),y1(j),rep,3),resp_var(x1(j),y1(j),rep,4),'Color',colors(idx_s(1,j),:),'Marker','.','MarkerSize',10);
    hold on
    if resp_var(x1(j),y1(j),rep,1)==-1
       plot3(resp_var(x1(j),y1(j),rep,2),resp_var(x1(j),y1(j),rep,3),resp_var(x1(j),y1(j),rep,4),'bo');
    end
    grid on
    xlabel('max tail deflection')
    xlim([0,120])
    %xlim([0,100])
    ylabel('min tail deflection')
    ylim([-120,0])
    %ylim([-100,0])
    zlabel('bidirectionality index')
    zlim([0,1])
    %             colors=colormap(hot(8))
    h=colorbar;
    h.Limits=[0,1]
    h.Ticks=1/(2*(ngroups+3)):1/(ngroups+3):(ngroups+2)/(ngroups+3)+1/(2*(ngroups+3))
    h.TickLabels={'1','2','3','4','5','6','7','8'};
end

%saveas(gcf,'ISIAll clusters mon 1.pdf')

%do the same for side 2
figure(72)
colors=colormap(hot(ngroups+3));

for j=1:length(x2)
    plot3(resp_var(x2(j),y2(j),rep,2),resp_var(x2(j),y2(j),rep,3),resp_var(x2(j),y2(j),rep,4),'Color',colors(idx_s(2,j),:),'Marker','.','MarkerSize',10);
    hold on
    grid on
    xlabel('max tail deflection')
    xlim([0,120])
    %xlim([0,100])
    ylabel('min tail deflection')
    ylim([-120,0])
    %ylim([-100,0])
    zlabel('bidirectionality index')
    zlim([0,1])
    %             colors=colormap(hot(8))
    h=colorbar;
    h.Limits=[0,1]
    h.Ticks=1/(2*(ngroups+3)):1/(ngroups+3):(ngroups+2)/(ngroups+3)+1/(2*(ngroups+3))
    h.TickLabels={'1','2','3','4','5','6','7','8'};
end
%saveas(gcf,'ISIAll clusters mon 2.pdf')    

%do stats on the escape time
 [P,ANOVATAB,STATS] = kruskalwallis([X(idx==2);X(idx==4);X(idx==5);X(idx==3)],[ones(size(X(idx==2)));2*ones(size(X(idx==4)));3*ones(size(X(idx==5)));4*ones(size(X(idx==3)))])
multcompare(STATS)

% save('../Classify_Escape_All.mat')
% CN=input('Please enter the cluster number corresponding to [ne, bi,u1,u2,sm]:');
% save('../Classify_Escape_All_final.mat')

%do a mesh plot for comparing fin movement and tail flicks
