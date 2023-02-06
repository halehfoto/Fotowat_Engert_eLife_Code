%This is a program that pools all the behavioral data, agnostic about
%imaging data. "Data for DMDL_DLDM_CBDL" under behavioral analysis under pooled data
clearvars
% close all
path=uigetdir;
cd(path)
load('ave_spont_rate.mat')
behav_fname1=dir('2*_b.mat');
stim_fname1=dir('2*_s.mat')
per2p=1.06295;
vidrate=25;
framerate=1/per2p;%2p frame rate determined based on resolution
trials1=1:1:10
trials=40*trials1;
trials20=40*(1:1:20);

%2p time
allflicks40=[]'
i=1;
for k=1:length(behav_fname1)
    load(behav_fname1(k).name)
    load(stim_fname1(k).name)
    if ISI==40
        %form a matrix with the first five columns corresponding to response =
        %1 and no response=0 and the sixth column corresponding to the ISI
        %if fps variable not present calculate it
        nframes=fpp;
        time=1/framerate:1/framerate:nframes/framerate;
        stim_time_trial1=time(round(ISI/per2p));
        ttc1=-1*lov/tand(theta_i);

        exist fps;
        if ~ans
            fps=fpp/reps;
        end
        exist stim_type;
        if ~ans
            ST(i)=1;
        else
            if strmatch(stim_type{1},'DL')
                ST(i)=1;%DL
            elseif strmatch(stim_type{1},'CB')
                ST(i)=2;%CB
            else
                ST(i)=3;%DM
            end
        end
        %Resample stimulus trace to match the length of flicks. This needs to
        %be done only for the first trial
        stim_temp=stim_size_degrees(1:fps);
        timev=1/vidrate:1/vidrate:length(flick(1,:))/vidrate;
        stim=interp1(time(1:fps),stim_temp, timev);
        %take the first 5 trials
    %     Normfact=max(max(flick(1:5,:),[],2));
        flicktemp=[];
        for j=1:20
            if size(flick,1)>=j
                tail_temp=fillmissing(flick(j,:),'spline');
                %find out if there was a response while the stimulus was more than
                %10 degrees and up to 5 s after it stops expanding
                %find stimulus values 
                [resp_time, Pmax, Pmin, PmaxminR, resp_amp, pm_amp, resp_dir,tail_dir_s,resp_temp]= quant_behav(flick(j,:),time,stim_size_degrees,ISI);
                resp(i,j)=resp_temp;
                resp_A(i,j)=abs(resp_amp)/50;
                resp_time_relcol(i,j)=ttc1/1000+resp_time- stim_time_trial1;
                %stim_resp(i,j)=trapz(abs(detrend(tail)))/length(index);
        %         if abs(stim_resp(i,j))>3
        %             stim_resp_bin(i,j)=sign(stim_resp(i,j))
        %         else
        %             stim_resp_bin(i,j)=0;
        %         end
        %         flicktemp=[flicktemp,(tail_temp(index))];
            else
                resp(i,j)=NaN;
                resp_A(i,j)=NaN;
                resp_time_relcol(i,j)=NaN;
            end

        end
    %     flicktemp=zscore((flicktemp));
    %     flicktemp_res=reshape(flicktemp,[length(tail),5]);
    %     for j=1:10
    %         if ISI==60
    %             allflicks60=[allflicks60;flicktemp_res(:,j)'];
    %         elseif ISI==40
    %             allflicks40=[allflicks40;flicktemp_res(:,j)'];
    %         else
    %             allflicks20=[allflicks20;flicktemp_res(:,j)'];
    %         end
    %     end
    %     
    clear fps
    clear nframes
    clear reps
    clear flick
    clear stim_type

        i=i+1;
    end
    
end


%% plot response probabilities fit the 10 trials
%subtract spontaneous rate for l/v=480 which is 0.1538
temp1=mean(resp(ST==1,1:10))-(1-poisscdf(0,0.0241*15.4231));;
temp1(temp1<0)=0;
figure
% plot(trials,temp1,'ro')
% hold on
% f1=fit(trials',temp1','exp1');
% p11 = predint(f1,trials',0.95,'observation','off');
% y1=f1.a*exp(f1.b*trials);
% h=plot(f1,40*(1:1:20),f1.a*exp(f1.b*(40*(1:1:20))))
%text(100,0.5,strcat('tau_DL=',num2str(1/f1.b,'%.2f')))
%legend off
%hold on
temp2=mean(resp(ST==2,1:10))-(1-poisscdf(0,Ave_Spont_resp_rate_CBDL*15.4231));%0.0241 is the spontaneous flick rate calculated using Calculate_First_trial_resp, and 15.4231 is the duration of the window used to detect flicks ;
temp2(temp2<0)=0;

temp3=mean(resp(ST==2,11:20))-(1-poisscdf(0,Ave_Spont_resp_rate_CBDL*15.4231));
temp3(temp3<0)=0;

plot(trials,temp2,'mo-')
hold on
plot(trials20(11:20)',temp3,'ro-')
f1=fit(trials',temp2','exp1');
f2=fit(trials20(11:20)',temp3','exp1');
y1=f1.a*exp(f1.b*trials);
y2=f2.a*exp(f2.b*trials20(11:20));

% h1=plot(f1,40*(1:1:10),f1.a*exp(f1.b*(40*(1:1:10))))
% h2=plot(f2,40*(11:1:20),f2.a*exp(f2.b*(40*(11:1:20))))

% plot(trials',p11,'k:')
% set(h,'color','r')
text(100,0.3,strcat('tau_CB=',num2str(1/f1.b,'%.2f')))
legend off
hold on

text(500,0.3,strcat('tau_DL=',num2str(1/f2.b,'%.2f')))
legend off
hold on
% set(h1,'color','m')
% set(h2,'color','r')
% plot(trials,mean(resp(ST==2,1:10)),'mo-')

% figure
% plot(trials(1:5),mean(resp(ST==3,1:5)),'ko-')
% hold on
% plot(trials(6:10),mean(resp(ST==3,6:10)),'ro-')
% ylim([0,0.9])

for i=1:10
    temp=sum(resp(ST==2,i));
    sum_all_flicks=round(temp);
    [PHAT, PCI] = binofit(sum_all_flicks,length(resp(ST==2,i)),0.05)
    pci_cb(1,i)=PCI(1);
    pci_cb(2,i)=PCI(2);

    phat_cb(i)=PHAT;
end
phat_cb_bs=phat_cb-(1-poisscdf(0,Ave_Spont_resp_rate_CBDL*15.4231));
phat_cb_bs(phat_cb_bs<0)=0;

plot(trials,phat_cb_bs,'mo-')
hold on

pci_cb_bs(1,:)=pci_cb(1,:)-(1-poisscdf(0,Ave_Spont_resp_rate_CBDL*15.4231));
pci_cb_bs(2,:)=pci_cb(2,:)-(1-poisscdf(0,Ave_Spont_resp_rate_CBDL*15.4231));

pci_cb_bs(1,pci_cb_bs(1,:)<0)=0;
pci_cb_bs(2,pci_cb_bs(2,:)<0)=0;

plot(trials,pci_cb_bs(1,:),'m:')
plot(trials,pci_cb_bs(2,:),'m:')


for i=11:20
    [PHAT, PCI] = binofit(sum(resp(ST==2,i)),length(resp(ST==2,i)),0.05)
    pci_dl2(1,i-10)=PCI(1);
    pci_dl2(2,i-10)=PCI(2);

    phat_dl2(i-10)=PHAT;
end

phat_dl2_bs=phat_dl2-(1-poisscdf(0,Ave_Spont_resp_rate_CBDL*15.4231));
phat_dl2_bs(phat_dl2_bs<0)=0;

pci_dl2_bs(1,:)=pci_dl2(1,:)-(1-poisscdf(0,Ave_Spont_resp_rate_CBDL*15.4231));
pci_dl2_bs(2,:)=pci_dl2(2,:)-(1-poisscdf(0,Ave_Spont_resp_rate_CBDL*15.4231));

pci_dl2_bs(1,pci_dl2_bs(1,:)<0)=0;
pci_dl2_bs(2,pci_dl2_bs(2,:)<0)=0;


plot(trials20(11:20)',phat_dl2_bs,'ro-')
hold on
plot(trials20(11:20)',pci_dl2_bs(1,:),'r:')
plot(trials20(11:20),pci_dl2_bs(2,:),'r:')
axis tight
box off

%% plot response probabilities all in one figure, latest version for the revision
figure(20)
hold on
subplot(2,1,2)
plot(trials20',[temp2,temp3],'mo-','LineWidth',1)
hold on
plot(trials20(11:20)',phat_dl2_bs,'ro-')
plot(trials20(11:20)',pci_dl2_bs(1,:),'r:')
plot(trials20(11:20),pci_dl2_bs(2,:),'r:')

% x2=[trials20(11:20),fliplr(trials20(11:20))];
% inBetween=[pci_dl2_bs(1,:), fliplr(pci_dl2_bs(2,:))];
% fill(x2,inBetween,'r','FaceAlpha',0.05)
% 
% x2=[trials,fliplr(trials)];
% inBetween=[pci_cb_bs(1,:), fliplr(pci_cb_bs(2,:))];
% fill(x2,inBetween,'m','FaceAlpha',0.05)
% 
plot(trials,pci_cb_bs(1,:),'m:')
plot(trials,pci_cb_bs(2,:),'m:')

x2=[trials20(11:20),fliplr(trials20(11:20))];
inBetween=[pci_dl2_bs(1,:), fliplr(pci_dl2_bs(2,:))];
fill(x2,inBetween,'r','FaceAlpha',0.05,'LineStyle','none')

x2=[trials20,fliplr(trials20)];
inBetween=[[pci_cb_bs(1,:),pci_dl2_bs(1,:)], fliplr([pci_cb_bs(2,:),pci_dl2_bs(2,:)])];
fill(x2,inBetween,'m','FaceAlpha',0.05,'LineStyle','none')
xlim([40,800])
ylim([-0.01,0.48])
box off

%% do exactly the same thing as above for response amplitude

figure
plot(trials,mean(resp_A(ST==1,1:10)),'ro-')
hold on
errorbar(trials,mean(resp_A(ST==1,1:10)),std(resp_A(ST==1,1:10))/sqrt(size(resp_A(ST==1),2)))
hold on
f1=fit(trials',(mean(resp_A(ST==1,1:10)))','exp1');
p11 = predint(f1,trials',0.95,'observation','off');
y1=f1.a*exp(f1.b*trials);
h=plot(f1,40*(1:1:10),f1.a*exp(f1.b*(40*(1:1:10))))
% plot(trials',p11,'k:')
set(h,'color','k')
text(100,1,strcat('tau_DL=',num2str(1/f1.b,'%.2f')))
legend off
hold on

h=plot(trials,mean(resp_A(ST==2,1:10)),'mo-')
errorbar(trials,mean(resp_A(ST==2,1:10)),std(resp_A(ST==2,1:10))/sqrt(size(resp_A(ST==2),2)),'Color','m')

hold on
plot(trials20(11:20)',mean(resp_A(ST==2,11:20)),'ro-')
errorbar(trials20(11:20),mean(resp_A(ST==2,11:20)),std(resp_A(ST==2,11:20))/sqrt(size(resp_A(ST==2),2)),'Color','r')

f1=fit(trials',(mean(resp_A(ST==2,1:10)))','exp1');
f2=fit(trials20(11:20)',(mean(resp_A(ST==2,11:20)))','exp1');
y1=f1.a*exp(f1.b*trials);
y2=f2.a*exp(f2.b*trials20(11:20));

h1=plot(f1,40*(1:1:10),f1.a*exp(f1.b*(40*(1:1:10))))
set(h,'color','k')

h2=plot(f2,40*(11:1:20),f2.a*exp(f2.b*(40*(11:1:20))))
set(h,'color','k')


% plot(trials',p11,'k:')
% set(h,'color','r')
text(100,0.4,strcat('tau_CB=',num2str(1/f1.b,'%.2f')))
legend off
hold on

text(500,0.4,strcat('tau_DL=',num2str(1/f2.b,'%.2f')))
legend off
hold on
% plot(trials,mean(resp_A(ST==2,1:10)),'mo-')




%% do the same for dim first
figure

plot(trials,mean(resp_A(ST==1,1:10)),'ro-')
hold on
errorbar(trials,mean(resp_A(ST==1,1:10)),std(resp_A(ST==1,1:10))/sqrt(size(resp_A(ST==1),2)))
hold on
f1=fit(trials',(mean(resp_A(ST==1,1:10)))','exp1');
p11 = predint(f1,trials',0.95,'observation','off');
y1=f1.a*exp(f1.b*trials);
h=plot(f1,40*(1:1:10),f1.a*exp(f1.b*(40*(1:1:10))))
% plot(trials',p11,'k:')
set(h,'color','k')
text(100,1,strcat('tau_DL=',num2str(1/f1.b,'%.2f')))
legend off
hold on
plot(trials(1:5),mean(resp_A(ST==3,1:5)),'ko-')
hold on
errorbar(trials(1:5),mean(resp_A(ST==3,1:5)),std(resp_A(ST==3,1:5))/sqrt(size(resp_A(ST==3),2)),'Color','k')

hold on
plot(trials(6:10),mean(resp_A(ST==3,6:10)),'ro-')
hold on
errorbar(trials(6:10),mean(resp_A(ST==3,6:10)),std(resp_A(ST==3,6:10))/sqrt(size(resp_A(ST==3),2)),'Color','r')
errorbar(trials(6:10),mean(resp_A(ST==3,6:10)),std(resp_A(ST==3,6:10)),'Color','r')
ylim([0,1.8])

%probability
figure
temp1=mean(resp(ST==1,1:10))-(1-poisscdf(0,0.0241*15.4231));
temp1(temp1<0)=0;

plot(trials,temp1,'ro-')
hold on
f1=fit(trials',temp1','exp1');
p11 = predint(f1,trials',0.95,'observation','off');
y1=f1.a*exp(f1.b*trials);
h=plot(f1,40*(1:1:10),f1.a*exp(f1.b*(40*(1:1:10))))
% plot(trials',p11,'k:')
% set(h,'color','r')
text(100,0.4,strcat('tau_DL=',num2str(1/f1.b,'%.2f')))
legend off
hold on

% plot(trials(1:5),mean(resp(ST==3,1:5)),'ko-')
% hold on
% plot(trials(6:10),mean(resp(ST==3,6:10)),'ro-')
% ylim([0,0.9])
%plot the corresponding confidence intervals
for i=1:5
    [PHAT, PCI] = binofit(sum(resp(ST==3,i)),length(resp(ST==3,i)),0.05)
    pci_dim(1,i)=PCI(1);
    pci_dim(2,i)=PCI(2);

    phat_dim(i)=PHAT;
end

phat_dim_bs=phat_dim-(1-poisscdf(0,Ave_Spont_resp_rate_DMDL*15.4231));
phat_dim_bs(phat_dim_bs<0)=0;
plot(trials(1:5),phat_dim_bs,'ko-')
hold on
pci_dim_bs(1,:)=pci_dim(1,:)-(1-poisscdf(0,Ave_Spont_resp_rate_DMDL*15.4231));
pci_dim_bs(1,pci_dim_bs(1,:)<0)=0;
pci_dim_bs(2,:)=pci_dim(2,:)-(1-poisscdf(0,Ave_Spont_resp_rate_DMDL*15.4231));
pci_dim_bs(2,pci_dim_bs(2,:)<0)=0;

plot(trials(1:5),pci_dim_bs(1,:),'k:')
plot(trials(1:5),pci_dim_bs(2,:),'k:')
for i=6:10
    [PHAT, PCI] = binofit(sum(resp(ST==3,i)),length(resp(ST==3,i)),0.05)
    pci_dl(1,i-5)=PCI(1);
    pci_dl(2,i-5)=PCI(2);

    phat_dl(i-5)=PHAT;
end
phat_dl_bs=phat_dl-(1-poisscdf(0,Ave_Spont_resp_rate_DMDL*15.4231));
phat_dl_bs(phat_dl_bs<0)=0;

pci_dl_bs(1,:)=pci_dl(1,:)-(1-poisscdf(0,Ave_Spont_resp_rate_DMDL*15.4231));
pci_dl_bs(1,pci_dl_bs(1,:)<0)=0;
pci_dl_bs(2,:)=pci_dl(2,:)-(1-poisscdf(0,Ave_Spont_resp_rate_DMDL*15.4231));
pci_dl_bs(2,pci_dl_bs(2,:)<0)=0;

plot(trials(6:10),phat_dl_bs,'r*-')
hold on
plot(trials(6:10),pci_dl_bs(1,:),'r:')
plot(trials(6:10),pci_dl_bs(2,:),'r:')
axis tight
box off
%% 
figure(20)
hold on
subplot(2,1,1)
plot(trials20(1:10)',[phat_dim_bs,phat_dl_bs],'bo-','LineWidth',1)
hold on
plot(trials20(6:10),phat_dl_bs,'ro-','LineWidth',1)
hold on
plot(trials20(1:5),pci_dim_bs(1,:),'b:')
plot(trials20(1:5),pci_dim_bs(2,:),'b:')
% axis tight
% ylim([-0.05,0.45])

x2=[trials20(1:10),fliplr(trials20(1:10))];
inBetween=[[pci_dim_bs(1,:),pci_dl_bs(1,:)], fliplr([pci_dim_bs(2,:),pci_dl_bs(2,:)])];
fill(x2,inBetween,'k','FaceAlpha',0.05,'LineStyle','none')

x2=[trials20(6:10),fliplr(trials20(6:10))];
inBetween=[pci_dl_bs(1,:), fliplr(pci_dl_bs(2,:))];
fill(x2,inBetween,'r','FaceAlpha',0.05,'LineStyle','none')
ylim([-0.01,0.48])
xlim([40,800])
box off
% hold on
plot(trials20(6:10),pci_dl_bs(1,:),'r:')
plot(trials20(6:10),pci_dl_bs(2,:),'r:')
% axis tight
% ylim([-0.05,0.45])
% x2=[trials20(10:15),fliplr(trials20(10:15))];
% inBetween=[pci_dl_bs(1,:), fliplr(pci_dl_bs(2,:))];
% fill(x2,inBetween,'r','FaceAlpha',0.05)
% ylim([0,0.45])
% box off
% 
%% generate bar plot to compare the model with the data from experiments from other fish
f1=fit(trials',(mean(resp(ST==1,1:10)))','exp1');
p11 = predint(f1,trials',0.95,'observation','off');
t=40*(1:1:20)
y1=f1.a*exp(f1.b*t);



X=categorical({'DL1','DL6','5XDM>DL'})
labels=reordercats(X,{'DL1','DL6','5XDM>DL'})
resp_prob_ave=[y1(1),y1(6),mean(resp(ST==3,6))]

%% generate bar plot to compare the model with the data from experiments from other fish
f1=fit(trials',(mean(resp(ST==1,1:10)))','exp1');
p11 = predint(f1,trials',0.95,'observation','off');
t=40*(1:1:20)
y1=f1.a*exp(f1.b*t);



X=categorical({'DL1','DL11','10xCB>DL'})
labels=reordercats(X,{'DL1','DL11','10xCB>DL'})
resp_prob_ave=[y1(1),y1(11),mean(resp(ST==2,11))]

figure
bar(labels,resp_prob_ave)
ylabel('resp_prob')
%%

load('DMDLCB_CrossHab_rates.mat');

plot(trials(6),mean(resp_all(Stim_type_all==2,2)),'ko');
plot(trials(10)+40,mean(resp_all(Stim_type_all==3,3)),'mo');
[phat,pci]=binofit(sum(resp_all(Stim_type_all==2,2)),length(resp_all(Stim_type_all==2,2)),0.05)
errorbar(trials(6),mean(resp_all(Stim_type_all==2,2)),pci(1)-mean(resp_all(Stim_type_all==2,2)),pci(2)-mean(resp_all(Stim_type_all==2,2)));

[phat,pci]=binofit(sum(resp_all(Stim_type_all==3,3)),length(resp_all(Stim_type_all==3,3)),0.05)
errorbar(trials(10)+40,mean(resp_all(Stim_type_all==3,3)),pci(1)-mean(resp_all(Stim_type_all==3,3)),pci(2)-mean(resp_all(Stim_type_all==3,3)));

nDLF=sum(Stim_type_all==1)
nDMF=sum(Stim_type_all==2)
nCBF=sum(Stim_type_all==3)


