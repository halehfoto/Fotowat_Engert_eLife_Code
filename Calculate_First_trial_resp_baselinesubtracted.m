%this is a program that calculates response to the first encounter with the stimulus
%for DL, DM and CB. Choose DLDMCB_trial1_comp. Created by Haleh Fotowat
clearvars
close all
path=uigetdir;
cd(path)
filename_b1=dir('2*_b.mat')
filename_s1=dir('2*_s.mat')
%DL: type 1, DM: type 2, CB: type 3
Stim_type_all=zeros(1,length(filename_s1));
per2p=1.06295;
framerate=1/per2p;%2p frame rate determined based on resolution
resp_all_trials_ISI=[];
resp_all_trials=[];
resp_time_all_trials_ISI=[];
all_spont_resp_dl=[];
all_spont_time_dl=[];
all_evoked_resp_dl=[];
all_evoked_time_dl=[];
all_spont_respcnt_dl=[];
k=1;
for i=1:length(filename_b1)
    filename_s1(i).name
    load(filename_s1(i).name);
    load(filename_b1(i).name);
    %     stim_type
    if ~exist('repse','var')
        repse=20;
    end
    if ~exist('stim_type','var')
        stim_type{1}='DL';
    end
    if strcmp(stim_type{1},'DL')
        Stim_type_all(i)=1;
    elseif strcmp(stim_type{1},'DM')
        Stim_type_all(i)=2;
    elseif strcmp(stim_type{1},'CB')
        Stim_type_all(i)=3;
    end
    nframes=fpp;
    fps=fpp/reps;
    ttc1=-1*lov/tand(theta_i);
    TTC(i)=ttc1;
    time=1/framerate:1/framerate:nframes/framerate;
    for j=1:1
        [resp_time, Pmax, Pmin, PmaxminR, resp_amp, pm_amp,resp_dir,tail_dir_s, resp, pm_amp_isi,resp_dir_isi, resp_isi, resp_amp_isi, resp_time_isi,time_window]=quant_behav(flick(j,:),time,stim_size_degrees,ISI);
        stim_time_trial=time(round(ISI/per2p):fps);
        STT(i)=stim_time_trial(1);
        resp_time_relstim(i,j)=resp_time-stim_time_trial(1);
        resp_time_relstart(i,j)=resp_time-time(find(stim_size_degrees>5,1));
        resp_time_relcol(i,j)=ttc1/1000+resp_time_relstim(i,j);
        abs_resp_amp_all(i,j)=abs(pm_amp(1));%amplitude of the first peak
        resp_amp_all(i,j)=(pm_amp(1));
        resp_dir_all(i,j)=resp_dir;%direction of the first peak
        resp_all(i,j)=resp;
        resp_time_all(i,j)=resp_time;
        
        % do the same for isi time_window
        TW(i)=time_window;
        %         resp_time_ISI(i,j)=resp_time_isi-stim_time_trial(1);%this is just relative to the first data point
        resp_time_ISI(i,j)=resp_time_isi-ISI+TW(i);
        resp_time_relcol_ISI(i,j)=ttc1/1000+resp_time_ISI(i,j);
        abs_resp_amp_all_isi(i,j)=abs(pm_amp_isi(1));%amplitude of the first peak
        resp_amp_all_isi(i,j)=(pm_amp_isi(1));
        resp_dir_all_isi(i,j)=resp_dir_isi;%direction of the first peak
        resp_all_isi(i,j)=resp_isi;
        
        
        if isnan(resp_time)
            rtime=-1;
        else
            rtime=resp_time_relstim(i,j);
        end
        behav_clust_param(i,:,j)=[rtime, Pmax, Pmin,PmaxminR];
        all_flicks(i,j,:)=tail_dir_s(25*ISI+1:end);
    end
    if Stim_type_all(i)==1
        if size(flick,1)>10 &ISI==40
            for j=1:10
                [resp_time, Pmax, Pmin, PmaxminR, resp_amp, pm_amp,resp_dir,tail_dir_s, resp, pm_amp_isi,resp_dir_isi, resp_isi, resp_amp_isi, resp_time_isi,resp_ISI_cnt,time_window]=quant_behav_ISIrespcnt(flick(j,:),time,stim_size_degrees,ISI);
                resp_all_trials_ISI=[resp_all_trials_ISI;resp_isi];
                resp_time_all_trials_ISI=[resp_time_all_trials_ISI, resp_time_isi-ISI+TW(i)];
                resp_per_trial_ISI(k,j)=resp_isi;
                resp_rate_per_ISI(k,j)=resp_ISI_cnt/TW(i);
                resp_per_trial(k,j)=resp;
                resp_per_time_trial_ISI(k,j)=resp_time_isi-ISI+TW(i);
                resp_per_time_trial(k,j)=resp_time-stim_time_trial(1)
%                 figure;
%                 plot(tail_dir_s);
%                 resp_ISI_cnt

            end
            k=k+1;
        end
    end
    if Stim_type_all(i)==1
        for kk=1:repse
            [resp_time, Pmax, Pmin, PmaxminR, resp_amp, pm_amp,resp_dir,tail_dir_s, resp, pm_amp_isi,resp_dir_isi, resp_isi, resp_amp_isi, resp_time_isi,resp_ISI_cnt,time_window]=quant_behav_ISIrespcnt(flick(j,:),time,stim_size_degrees,ISI);
            all_spont_resp_dl=[all_spont_resp_dl,resp_isi];
            all_spont_respcnt_dl=[all_spont_respcnt_dl,resp_ISI_cnt/TW(i)];
            all_spont_time_dl=[all_spont_time_dl,resp_time_isi-ISI+TW(i)];
            all_evoked_resp_dl=[all_evoked_resp_dl,resp];
            all_evoked_time_dl=[all_evoked_time_dl,resp_time-stim_time_trial(1)];

        end
    end
        
%     resp_amp_all(i,1)
%     resp_dir_all(i,1)
%     (resp_time_relstim(i,1)+stim_time_trial(1))*25
    clear nframes
    clear fps
    clear reps
    clear flick
    clear stim_type
    %close all
    
end

figure;
hold on; 
subplot(2,1,1)
plot(nanmean(resp_per_trial_ISI),'bo-');
hold on
plot(mean(resp_per_trial),'ro-')
Ave_Spont_resp=mean(nanmean(resp_per_trial_ISI))

ylabel('response probability')
xlabel('trial number')
subplot(2,1,2);stdshade(resp_rate_per_ISI,0.1,'k')

ylabel('spontaneous flick rate (Hz)')
xlabel('trial number')

figure;plot(nanmean(resp_per_time_trial_ISI),'bo-')
hold on
plot(nanmean(resp_per_time_trial),'ro-')
hold on

Ave_Spont_resp_rate=mean(nanmean(resp_rate_per_ISI))



[Ns,EDGES,BIN]=histcounts(reshape(resp_per_time_trial_ISI,1,numel(resp_per_time_trial_ISI)),7)
ylim([0,30])
figure;bar(Ns)
figure;hist(reshape(resp_per_time_trial_ISI,1,numel(resp_per_time_trial_ISI)),7)
figure;hist(reshape(resp_per_time_trial,1,numel(resp_per_time_trial)),7)

[Nl,EDGES,BIN]=histcounts(reshape(resp_per_time_trial,1,numel(resp_per_time_trial)),7)
figure;hist(reshape(resp_per_time_trial,1,numel(resp_per_time_trial)),7)
figure;bar(Nl)
figure;bar(Nl-nanmean(Ns))

figure;histogram(~isnan(reshape(resp_per_time_trial_ISI,1,numel(resp_per_time_trial_ISI))))



% resp_all_DL=[];
% for i=1:length(filename_b)    
%     load(filename_s(i).name)
%     load(filename_b(i).name)
%     if ~exist('stim_type','var')
%         stim_type{1}='DL';
%     end
%     if strcmp(stim_type{1},'DL')
%         Stim_type_all_DL(k)=1;
%         nframes=fpp;
%         fps=fpp/reps;
%         ttc1=-1*lov/tand(theta_i);
%         TTC(i)=ttc1;
%         time=1/framerate:1/framerate:nframes/framerate;
%         for j=1:5
%             [resp_time, Pmax, Pmin, PmaxminR, resp_amp, pm_amp,resp_dir,tail_dir_s, resp, pm_amp_isi,resp_dir_isi, resp_isi, resp_amp_isi, resp_time_isi,time_window]=quant_behav(flick(j,:),time,stim_size_degrees,ISI);
%             resp_all_DL(k,j)=resp;
%             
%         end
%         k=k+1
%     end
%     clear stim_type
% end
% %what is the probabilty of getting responses in the subsequent trials if
% %the first trial is zero
% index=find(resp_all_DL(:,1)==0);
% index2=find(resp_all_DL(:,1)==1);
% 
% index3=find(resp_all_DL(:,1)==0 & resp_all_DL(:,2)==0 );
% index4=find(resp_all_DL(:,1)==0 & resp_all_DL(:,2)==0 & resp_all_DL(:,3)==0);
% 
% for j=1:5
%     resp_prob_1_0(j)=sum(resp_all_DL(index,j))/length(index);
%     resp_prob_1_1(j)=sum(resp_all_DL(index2,j))/length(index2);
%     resp_prob_1_0_0(j)=sum(resp_all_DL(index3,j))/length(index3);
%     resp_prob_1_0_0_0(j)=sum(resp_all_DL(index4,j))/length(index4);
% 
% end
% index1=1:1:length(resp_all_DL);
% index1(index4)=[];
% sum(resp_all_DL(:,1))/length(index1)

X3=categorical({'DL','CB','DM','SPONT'})
labels3=reordercats(X3,{'DL','CB','DM','SPONT'})

figure;resp_prob=[(sum(resp_all(Stim_type_all'==1))/sum(Stim_type_all'==1))-(1-poisscdf(0,0.0241*15.4231)),(sum(resp_all(Stim_type_all'==3))/sum(Stim_type_all'==3))-(1-poisscdf(0,0.0241*15.4231)),(sum(resp_all(Stim_type_all'==2))/sum(Stim_type_all'==2))-(1-poisscdf(0,0.0241*15.4231)),(sum(~isnan(resp_dir_all_isi))/length(resp_dir_all_isi))-(1-poisscdf(0,0.0241*15.4231))]
bar(labels3,resp_prob)
hold on
title('Response Probability')

X3=categorical({'DL','CB','DM'})
labels3=reordercats(X3,{'DL','CB','DM'})

figure;resp_prob=[(sum(resp_all(Stim_type_all'==1))/sum(Stim_type_all'==1))-(1-poisscdf(0,0.0241*15.4231)),(sum(resp_all(Stim_type_all'==3))/sum(Stim_type_all'==3))-(1-poisscdf(0,0.0241*15.4231)),(sum(resp_all(Stim_type_all'==2))/sum(Stim_type_all'==2))-(1-poisscdf(0,0.0241*15.4231))]
bar(labels3,resp_prob)
hold on
title('Response Probability')



resp_time=[(resp_time_relcol(Stim_type_all==1));(resp_time_relcol(Stim_type_all==3));(resp_time_relcol(Stim_type_all==2));resp_time_relcol_ISI]
groups=[ones(size((resp_time_relcol(Stim_type_all==1))));2*ones(size((resp_time_relcol(Stim_type_all==3))));3*ones(size((resp_time_relcol(Stim_type_all==2))));4*ones(size(resp_time_relcol_ISI))]
[a,b,STATS]=kruskalwallis(resp_time,groups)
figure;multcompare(STATS)

resp_time=[(resp_time_relcol(Stim_type_all==1));(resp_time_relcol(Stim_type_all==3));(resp_time_relcol(Stim_type_all==2))]
groups=[ones(size((resp_time_relcol(Stim_type_all==1))));2*ones(size((resp_time_relcol(Stim_type_all==3))));3*ones(size((resp_time_relcol(Stim_type_all==2))))]
[a,b,STATS]=kruskalwallis(resp_time,groups)
figure;A=multcompare(STATS)


resp_time_relstart=[(resp_time_relstim(Stim_type_all==1));(resp_time_relstim(Stim_type_all==3));(resp_time_relstim(Stim_type_all==2));resp_time_ISI]
groups=[ones(size((resp_time_relcol(Stim_type_all==1))));2*ones(size((resp_time_relcol(Stim_type_all==3))));3*ones(size((resp_time_relcol(Stim_type_all==2))));4*ones(size(resp_time_ISI))]
[a,b,STATS]=kruskalwallis(resp_time_relstart,groups)
figure;A=multcompare(STATS)


%compare response amplitudes
resp_amp=[abs_resp_amp_all(Stim_type_all==1)/50;abs_resp_amp_all(Stim_type_all==3)/50;abs_resp_amp_all(Stim_type_all==2)/50;abs_resp_amp_all_isi/50];
groups=[ones(size(abs_resp_amp_all(Stim_type_all==1)));2*ones(size(abs_resp_amp_all(Stim_type_all==3)));3*ones(size(abs_resp_amp_all(Stim_type_all==2)));4*ones(size(abs_resp_amp_all_isi))]
[a,b,STATS]=kruskalwallis(resp_amp,groups)

resp_amp=[abs_resp_amp_all(Stim_type_all==1)/50;abs_resp_amp_all(Stim_type_all==3)/50;abs_resp_amp_all(Stim_type_all==2)/50];
groups=[ones(size(abs_resp_amp_all(Stim_type_all==1)));2*ones(size(abs_resp_amp_all(Stim_type_all==3)));3*ones(size(abs_resp_amp_all(Stim_type_all==2)))]
[a,b,STATS]=kruskalwallis(resp_amp,groups)
figure;A=multcompare(STATS)

resp_amp=[abs_resp_amp_all(Stim_type_all==1)/50;abs_resp_amp_all_isi/50];
groups=[ones(size(abs_resp_amp_all(Stim_type_all==1)));2*ones(size(abs_resp_amp_all_isi))]
[a,b,STATS]=kruskalwallis(resp_amp,groups)



resp_amp=[abs_resp_amp_all(Stim_type_all==1)/50;abs_resp_amp_all(Stim_type_all==3)/50;abs_resp_amp_all(Stim_type_all==2)/50];
groups=[ones(size(abs_resp_amp_all(Stim_type_all==1)));2*ones(size(abs_resp_amp_all(Stim_type_all==3)));3*ones(size(abs_resp_amp_all(Stim_type_all==2)))]
[a,b,STATS]=kruskalwallis(resp_amp,groups)
multcompare(STATS)
%save('First_trial_Resp.mat');
save('ave_spont_rate.mat','Ave_Spont_resp','Ave_Spont_resp_rate')