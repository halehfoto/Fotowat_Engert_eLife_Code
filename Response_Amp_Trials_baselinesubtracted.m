%This is a code to generate response amplitude changes with trial and side
%plots. 
close all
clearvars
s=what;
size(s.mat)
cnt=1;
for i=1:length(s.mat)
    filename=s.mat{i}
    L=length(filename);
    X=filename(L-9:L-8);
    
    if strmatch(X,'03')
        ISI_fish(cnt:cnt+9)=180
        cnt=cnt+10;
    elseif strmatch(X,'05')
        ISI_fish(cnt:cnt+9)=5
        cnt=cnt+10;
    else
        ISI_fish(cnt:cnt+9)=10
        cnt=cnt+10;
    end
end
load('../Classify_Escape_All_final.mat')
ne=CN(1);
bi=CN(2);
u1=CN(3);
u2=CN(4);
sm=CN(5);
idx_side=[idx(1:length(idx)/2),idx((length(idx)/2)+1:end)]
rep=1;
for i=1:length(s.mat)
    for trial=1:10
        for k=1:2
            resp_time_temp=resp_var((i-1)*10+trial,k,rep,10);
            if resp_time_temp*1000< ttc_stim_save{1}(find(ang_size_save{1}>5,1))|| resp_time_temp*1000>(ttc_stim_save{1}(find(ang_size_save{1}==max(ang_size_save{1}),1))+5000)
                flag((i-1)*10+trial,k,rep)=1;
                resp_var((i-1)*10+trial,k,rep,9)=0; %set the response to zero if it occurs outside of the correct window
                resp_var((i-1)*10+trial,k,rep,10)=NaN;
                resp_var((i-1)*10+trial,k,rep,11)=NaN;
 
            else
                flag((i-1)*10+trial,k,rep)=1;
            end
        end
    end
end

% k=1;
% rep=1;
close all
for k=1:2
    for trial=1:10
        index=find(ISI_fish'==180 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0);
        Ave_time_180(trial,k)=nanmean((resp_var(index,k,1,10)));
        Ste_time_180(trial,k)=nanstd((resp_var(index,k,1,10)))/sqrt(sum(~isnan(resp_var(index,k,1,10))));
        Ave_amp_180(trial,k)=nanmean(abs(resp_var(index,k,1,9)));
        Ste_amp_180(trial,k)=nanstd(abs(resp_var(index,k,1,9)))/sqrt(sum(~isnan(resp_var(index,k,1,9))));
       % Npeaks(index,k)-9.2583*0.0241
        temp=sum(~isnan(resp_var(index,k,1,10)));
        
        sum_all_flicks=round(temp);
        [phat,pci]=binofit(sum_all_flicks,length(index),0.05) % subtract spontaneous rate
        resp_prob_180(trial,k)=phat;
        conf_180(trial,k,:)=pci;
        P_c_180(trial,k)=sum(resp_var(index,k,1,11)==1)/sum(~isnan(resp_var(index,k,1,11)));
 
        index_bi=find(idx_side(:,k)==bi & ISI_fish'==180 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0);
        Ave_time_180_bi(trial,k)=nanmean((resp_var(index_bi,k,1,10)));
        Ste_time_180_bi(trial,k)=nanstd((resp_var(index_bi,k,1,10)))/sqrt(sum(~isnan(resp_var(index_bi,k,1,10))));

        index_uni=find((idx_side(:,k)==u1 | idx_side(:,k)==u2)  & ISI_fish'==180 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0);
        Ave_time_180_uni(trial,k)=nanmean((resp_var(index_uni,k,1,10)));
        Ste_time_180_uni(trial,k)=nanstd((resp_var(index_uni,k,1,10)))/sqrt(sum(~isnan(resp_var(index_uni,k,1,10))));
       
        index_sm=find((idx_side(:,k)==sm)  & ISI_fish'==180 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0);
        Ave_time_180_sm(trial,k)=nanmean((resp_var(index_sm,k,1,10)));
        Ste_time_180_sm(trial,k)=nanstd((resp_var(index_sm,k,1,10)))/sqrt(sum(~isnan(resp_var(index_sm,k,1,10))));
        tcnt_180(trial,k)=sum(~isnan(resp_var(index,k,1,11)));
        tcnt_all_180(trial,k)=length(index)
    end
end



for k=1:2
    for trial=1:10
        index=find(ISI_fish'==10 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0);
        Ave_time_10(trial,k)=nanmean((resp_var(index,k,1,10)));
        Ste_time_10(trial,k)=nanstd((resp_var(index,k,1,10)))/sqrt(sum(~isnan(resp_var(index,k,1,10))));
        Std_time_10(trial,k)=nanstd((resp_var(index,k,1,10)));
        Ave_amp_10(trial,k)=nanmean(abs(resp_var(index,k,1,9)));
        Ste_amp_10(trial,k)=nanstd(abs(resp_var(index,k,1,9)))/sqrt(sum(~isnan(resp_var(index,k,1,9))));
        temp=sum(~isnan(resp_var(index,k,1,10)))%;-9.2583*0.0241*length(index);%subtract spontaneous spikes
        sum_all_flicks=round(temp);
        [phat,pci]=binofit(sum_all_flicks,length(index),0.05)
        resp_prob_10(trial,k)=phat;
        conf_10(trial,k,:)=pci;
        P_c_10(trial,k)=sum(resp_var(index,k,1,11)==1)/sum(~isnan(resp_var(index,k,1,11)));
        tcnt_10(trial,k)=sum(~isnan(resp_var(index,k,1,11)));
       tcnt_all_10(trial,k)=length(index)

    end
end


for k=1:2
    for trial=1:10
        index=find(ISI_fish'==5 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0);
        Ave_time_5(trial,k)=nanmean((resp_var(index,k,1,10)));
        Ste_time_5(trial,k)=nanstd((resp_var(index,k,1,10)))/sqrt(sum(~isnan(resp_var(index,k,1,10))));

        Ave_amp_5(trial,k)=nanmean(abs(resp_var(index,k,1,9)));
        Ste_amp_5(trial,k)=nanstd(abs(resp_var(index,k,1,9)))/sqrt(sum(~isnan(resp_var(index,k,1,9))));
        temp=sum(~isnan(resp_var(index,k,1,10)));%-9.2583*0.0241*length(index);
        sum_all_flicks=round(temp);
%         if sum_all_flicks<0
%             sum_all_flicks=0;
%         end
        [phat,pci]=binofit(sum_all_flicks,length(index),0.05)
        resp_prob_5(trial,k)=phat;
        conf_5(trial,k,:)=pci;
        P_c_5(trial,k)=sum(resp_var(index,k,1,11)==1)/sum(~isnan(resp_var(index,k,1,11)));
        tcnt_5(trial,k)=sum(~isnan(resp_var(index,k,1,11)));

    end
end


%% plot response probability and time
figure(1)
trials1=1:1:10;
trials2=11:1:20;
for i=1:2
     subplot(4,1,2)
    if i==1
        trials=5*trials1;
        color='r'
    else
        trials=5*trials2;
        color='b'
    end
    plot(trials,Ave_amp_5(:,i),'Color',color,'Marker','o')
    hold on
    errorbar(trials,Ave_amp_5(:,i),Ste_amp_5(:,i),'Color',color,'Marker','o');
    f1=fit(trials',Ave_amp_5(:,i),'exp1');
    p11 = predint(f1,trials',0.95,'functional','off');
    confint_amp_5(:,1:2,i)=p11;

    plot(f1,trials',Ave_amp_5(:,i)), hold on, plot(trials',p11,'m--')

%     y1=f1.a*exp(f1.b*trials);
%     h=plot(f1,trials,Ave_amp_5(:,i))
%     set(h,'color','k')

    text((i-1)*5+10*(i-1)*5+5,0.5,strcat('tau=',num2str(1/f1.b,'%.2f')))
    legend off
        set(gca,'Xlabel',[]) 

    resp_prob_5_bs(:,i)=resp_prob_5(:,i)-(1-poisscdf(0,0.0241*9.2583));%0.0241 is the spontaneous flick rate calculated using Calculate_First_trial_resp, and 9.2583 is the duration of the window used to detect flicks 
   resp_prob_5_bs(resp_prob_5_bs(:,i)<0,i)=0;
    axis tight
    box off

    subplot(4,1,1)
    plot(trials,resp_prob_5_bs(:,i),'Color',color,'Marker','o')
    hold on
    %plot the binary confidence intervals
    conf_5(:,i,1)=conf_5(:,i,1)-(1-poisscdf(0,0.0241*9.2583));
    conf_5(:,i,2)=conf_5(:,i,2)-(1-poisscdf(0,0.0241*9.2583));
    conf_5(conf_5(:,i,1)<0,i,1)=0;
    conf_5(conf_5(:,i,2)<0,i,2)=0;
   
    plot(trials,conf_5(:,i,1),'Color',color,'LineStyle',':')
    plot(trials,conf_5(:,i,2),'Color',color,'LineStyle',':')
    f1=fit(trials',resp_prob_5_bs(:,i),'exp1');
    %this is nonsimultaneous 95% confidence interval on individual elemtns
    %and the 'observation' type measures the uncertainy of predicting the
    %fitted curve plus random variations in the new observation
    %p11 = predint(f1,trials',0.95,'functional','off');
    %confint_5(:,1:2,i)=p11;

    plot(f1,trials',resp_prob_5_bs(:,i)), hold on, %plot(trials',p11,'m--')
        %daspect([0.1 0.01 1])
%     a=daspect;
%     daspect([0.5*a(1) a(2) a(3)])
%     y1=f1.a*exp(f1.b*trials);
%     h=plot(f1,trials,resp_prob_5(:,i))
%     set(h,'color','k')
    text((i-1)*10+5*(i-1)*10+5,0.5,strcat('tau=',num2str(1/f1.b,'%.2f')))
    legend off
    set(gca,'Xlabel',[]) 
    axis tight
    box off

    subplot(4,1,3)
    plot(trials,Ave_time_5(:,i),'Color',color,'Marker','o')
    hold on
    errorbar(trials,Ave_time_5(:,i),Ste_time_5(:,i),'Color',color,'Marker','o');
        set(gca,'Xlabel',[]) 
    axis tight
    box off

    subplot(4,1,4)
    plot(trials,P_c_5(:,i),'Color',color,'Marker','o')
    hold on
       set(gca,'Xlabel',[]) 
    axis tight
    box off
% 
    

% 
end
figure(2)
for i=1:2
    subplot(4,1,2)
    if i==1
        trials=10*trials1;
        color='r'
    else
        trials=10*trials2;
        color='b'
    end
    plot(trials,Ave_amp_10(:,i),'Color',color,'Marker','o')
    hold on
    errorbar(trials,Ave_amp_10(:,i),Ste_amp_10(:,i),'Color',color,'Marker','o');
    f1=fit(trials',Ave_amp_10(:,i),'exp1');
    p11 = predint(f1,trials',0.95,'functional','off');
    confint_amp_10(:,1:2,i)=p11;

    plot(f1,trials',Ave_amp_10(:,i)), hold on, plot(trials',p11,'m--')


    y1=f1.a*exp(f1.b*trials);
    h=plot(f1,trials,Ave_amp_10(:,i))
    set(h,'color','k')

    text((i-1)*10+10*(i-1)*10+10,0.5,strcat('tau=',num2str(1/f1.b,'%.2f')))
    legend off
   resp_prob_10_bs(:,i)=resp_prob_10(:,i)-(1-poisscdf(0,0.0241*9.2583));%0.0241 is the spontaneous flick rate calculated using Calculate_First_trial_resp, and 9.2583 is the duration of the window used to detect flicks 
   resp_prob_10_bs(resp_prob_10_bs(:,i)<0,i)=0;
    conf_10(:,i,1)=conf_10(:,i,1)-(1-poisscdf(0,0.0241*9.2583));
    conf_10(:,i,2)=conf_10(:,i,2)-(1-poisscdf(0,0.0241*9.2583));
    conf_10(conf_10(:,i,1)<0,i,1)=0;
    conf_10(conf_10(:,i,2)<0,i,2)=0;

   
;%-0.153; %baseline flick subtracted
    set(gca,'Xlabel',[]) 
    axis tight
    box off
    subplot(4,1,1)
    plot(trials,resp_prob_10_bs(:,i),'Color',color,'Marker','o')
    hold on
    %plot the binary confidence intervals
    plot(trials,conf_10(:,i,1),'Color',color,'LineStyle',':')
    plot(trials,conf_10(:,i,2),'Color',color,'LineStyle',':')

    f1=fit(trials',resp_prob_10_bs(:,i),'exp1');
    p11 = predint(f1,trials',0.95,'functional','off');
    confint_10(:,1:2,i)=p11;
    plot(f1,trials',resp_prob_10_bs(:,i)), hold on, %plot(trials',p11,'m--')

%     y1=f1.a*exp(f1.b*trials);
%     h=plot(f1,trials,resp_prob_10(:,i))
%     set(h,'color','k')
    text((i-1)*10+10*(i-1)*10+10,0.5,strcat('tau=',num2str(1/f1.b,'%.2f')))
    legend off
    set(gca,'Xlabel',[]) 
        axis tight
    box off

    subplot(4,1,3)
    plot(trials,Ave_time_10(:,i),'Color',color,'Marker','o')
    hold on
    errorbar(trials,Ave_time_10(:,i),Ste_time_10(:,i),'Color',color,'Marker','o');
    set(gca,'Xlabel',[]) 
    axis tight
    box off
    subplot(4,1,4)
    plot(trials,P_c_10(:,i),'Color',color,'Marker','o')
    hold on
    set(gca,'Xlabel',[]) 
    axis tight
    box off
    

end
figure(3)
for i=1:2
    subplot(4,1,2)
    if i==1
        trials=180*trials1;
        color='r';
    else
        trials=180*trials2;
        color='b';
    end
    plot(trials,Ave_amp_180(:,i),'Color',color,'Marker','o')
    hold on
    errorbar(trials,Ave_amp_180(:,i),Ste_amp_180(:,i),'Color',color,'Marker','o');
    f1=fit(trials',Ave_amp_180(:,i),'exp1');
    p11 = predint(f1,trials',0.95,'functional','off');
    confint_amp_180(:,1:2,i)=p11;

    plot(f1,trials',Ave_amp_180(:,i)), hold on, plot(trials',p11,'m--')

%     y1=f1.a*exp(f1.b*trials);
%     h=plot(f1,trials,Ave_amp_180(:,i))
%     set(h,'color','k')

    text((i-1)*180+10*(i-1)*180+180,0.5,strcat('tau=',num2str(1/f1.b,'%.2f')))
    legend off
     axis tight
    box off

   resp_prob_180_bs(:,i)=resp_prob_180(:,i)-(1-poisscdf(0,0.0241*9.2583));%0.0241 is the spontaneous flick rate calculated using Calculate_First_trial_resp, and 9.2583 is the duration of the window used to detect flicks 
   resp_prob_180_bs(resp_prob_180_bs(:,i)<0,i)=0;
    conf_180(:,i,1)=conf_180(:,i,1)-(1-poisscdf(0,0.0241*9.2583));
    conf_180(:,i,2)=conf_180(:,i,2)-(1-poisscdf(0,0.0241*9.2583));
    conf_180(conf_180(:,i,1)<0,i,1)=0;
    conf_180(conf_180(:,i,2)<0,i,2)=0;

    subplot(4,1,1)
    plot(trials,resp_prob_180_bs(:,i),'Color',color,'Marker','o')
    hold on
    %plot the binary confidence intervals
    plot(trials,conf_180(:,i,1),'Color',color,'LineStyle',':')
    plot(trials,conf_180(:,i,2),'Color',color,'LineStyle',':')
     %% could not fit an exponential had to manually run the program to plot the other two panels.
    f1=fit(trials',resp_prob_180_bs(:,i),'exp1');
    
    p11 = predint(f1,trials',0.95,'functional','off');
    confint_180(:,1:2,i)=p11;

    plot(f1,trials',resp_prob_180_bs(:,i)), hold on, %plot(trials',p11,'m--')
%     y1=f1.a*exp(f1.b*trials);
%     h=plot(f1,trials,resp_prob_180(:,i))
%     set(h,'color','k')
    text((i-1)*180+180*(i-1)*10+180,0.5,strcat('tau=',num2str(1/f1.b,'%.2f')))
    axis tight
    box off

    legend off
    subplot(4,1,3)
    plot(trials,Ave_time_180(:,i),'ro-')
    hold on
    errorbar(trials,Ave_time_180(:,i),Ste_time_180(:,i),'Color',color,'Marker','o');
    axis tight
    box off

    subplot(4,1,4)
    plot(trials,P_c_180(:,i),'Color',color,'Marker','o')
    hold on
    axis tight
    box off
end


%% plot response directions

figure
for k=1:2
    k1=1;
    m1=1;

    if k==1
        color='r';
    else
        color='b'
    end
    for trial=1:10
        k1=1;
        m1=1;

        index=find(ISI_fish'==180 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0);

        for j=1:length(index)
            if resp_var(index(j),k,1,11)==1
                plot((k-1)*10+trial,k1,'Color',color,'Marker','.','markersize',20)
                k1=k1+1;
                hold on
            elseif resp_var(index(j),k,1,11)==0
                plot((k-1)*10+trial,-1*m1,'Color',color,'Marker','o')
                m1=m1+1;
                hold on
            else
                plot((k-1)*10+trial,-1*m1,'Color',color,'Marker','x')
                m1=m1+1;
                hold on
                
            end
        end
        P_c_180(trial,k)=sum(resp_var(index,k,1,11)==1)/sum(~isnan(resp_var(index,k,1,11)));

    end
end

xlim([0,21])
ylim([-15,15])
line([0,21],[0,0])
xlabel('Trial number')
ylabel('Response direction')
title(strcat('ISI = 3 min,',num2str(sum((ISI_fish==180))/10),' fish'))
figure
for k=1:2
    k1=1;
    m1=1;

    if k==1
        color='r';
    else
        color='b'
    end
    for trial=1:10
        k1=1;
        m1=1;

        index=find(ISI_fish'==10 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0 );

        for j=1:length(index)
            if resp_var(index(j),k,1,11)==1
                plot((k-1)*10+trial,k1,'Color',color,'Marker','.','markersize',20)
                k1=k1+1;
                hold on
            elseif resp_var(index(j),k,1,11)==0
                plot((k-1)*10+trial,-1*m1,'Color',color,'Marker','o')
                m1=m1+1;
                hold on
            else
                plot((k-1)*10+trial,-1*m1,'Color',color,'Marker','x')
                m1=m1+1;
                hold on
                
            end
        end
        P_c_10(trial,k)=sum(resp_var(index,k,1,11)==1)/sum(~isnan(resp_var(index,k,1,11)));

    end
end
[a,b]=corrcoef(P_c_10(:,1),1:1:10)
[a,b]=corrcoef(P_c_10(:,2),1:1:10)
mean(P_c_10(:,1))
std(P_c_10(:,1))

xlim([0,21])
ylim([-20,20])
line([0,21],[0,0])
xlabel('Trial number')
ylabel('Response direction')
title(strcat('ISI = 10 s,',num2str(sum((ISI_fish==10))/10),' fish'))
figure
for k=1:2
    k1=1;
    m1=1;

    if k==1
        color='r';
    else
        color='b'
    end
    for trial=1:10
        k1=1;
        m1=1;

        index=find(ISI_fish'==5 & resp_var(:,k,1,6)==trial & flag(:,k,1)~=0);

        for j=1:length(index)
            if resp_var(index(j),k,1,11)==1
                plot((k-1)*10+trial,k1,'Color',color,'Marker','.','markersize',20)
                k1=k1+1;
                hold on
            elseif resp_var(index(j),k,1,11)==0
                plot((k-1)*10+trial,-1*m1,'Color',color,'Marker','o')
                m1=m1+1;
                hold on
            else
                plot((k-1)*10+trial,-1*m1,'Color',color,'Marker','x')
                m1=m1+1;
                hold on
                
            end
        end
        P_c_5(trial,k)=sum(resp_var(index,k,1,11)==1)/sum(~isnan(resp_var(index,k,1,11)));

    end
end
[a,b]=corrcoef(P_c_5(:,1),1:1:10)
[a,b]=corrcoef(P_c_5(:,2),1:1:10)
mean(P_c_5(:,1))
std(P_c_5(:,1))


xlim([0,21])
ylim([-15,15])
line([0,21],[0,0])
xlabel('Trial number')
ylabel('Response direction')
title(strcat('ISI =5 s,',num2str(sum((ISI_fish==5))/10),' fish'))
% 
% 
amps_180=[abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0,1,1,9));abs(resp_var(ISI_fish'==180 & flag(:,2,1)~=0,2,1,9))];
Group_180=[resp_var(ISI_fish'==180 & flag(:,1,1)~=0,1,1,6);10+resp_var(ISI_fish'==180 & flag(:,2,1)~=0,2,1,6)];
[P,ANOVATAB,STATS]=kruskalwallis(amps_180, Group_180);
figure
multcompare(STATS)


% %compare first three and last three trials
amps_180_3=[abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)<4 ,1,1,9));abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)>7,1,1,9));abs(resp_var(ISI_fish'==180 & flag(:,2,1)~=0& resp_var(:,2,1,6)<4,2,1,9))];
Group_180_3=[ones(size(abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)<4 ,1,1,9))));2*ones(size(abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)>7,1,1,9))));3*ones(size(abs(resp_var(ISI_fish'==180 & flag(:,2,1)~=0& resp_var(:,2,1,6)<4,2,1,9))))];
[P,ANOVATAB,STATS]=kruskalwallis(amps_180_3, Group_180_3);
figure
multcompare(STATS)
% %compare first three and last one and first on the other side

amps_180_3=[abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9));abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==2 ,1,1,9));abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==3 ,1,1,9));abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9));abs(resp_var(ISI_fish'==180 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9))];
Group_180_3=[ones(size(abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9))));2*ones(size(abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==2 ,1,1,9))));3*ones(size(abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==3 ,1,1,9))));4*ones(size(abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9))));5*ones(size(abs(resp_var(ISI_fish'==180 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9))))];;
[P,ANOVATAB,STATS]=kruskalwallis(amps_180_3, Group_180_3);
figure
multcompare(STATS)


%compare first, last and first

X1=abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9));
X10=abs(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9));
Y1=abs(resp_var(ISI_fish'==180 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9));
%first side one and first side 2
kruskalwallis([X1;Y1],[ones(size(X1));2*ones(size(Y1))])
kruskalwallis([X1;X10],[ones(size(X1));2*ones(size(X10))])
kruskalwallis([X10;Y1],[ones(size(X10));2*ones(size(Y1))])

%compare first, last and first

% amps_10_flf=[abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9));abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9));abs(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9))];
% Group_10_flf=[ones(size(abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9))));2*ones(size(abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9))));3*ones(size(abs(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9))))];
X1=abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9));
X2=abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==2 ,1,1,9));
X3=abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==3 ,1,1,9));
X10=abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9));
Y1=abs(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9));
%first side one and first side 2
kruskalwallis([X1;Y1],[ones(size(X1));2*ones(size(Y1))])
kruskalwallis([X2;Y1],[ones(size(X10));2*ones(size(Y1))])
kruskalwallis([X3;Y1],[ones(size(X10));2*ones(size(Y1))])

kruskalwallis([X1;X10],[ones(size(X1));2*ones(size(X10))])
kruskalwallis([X10;Y1],[ones(size(X10));2*ones(size(Y1))])


%for 5
X1=abs(resp_var(ISI_fish'==5 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9));
X2=abs(resp_var(ISI_fish'==5 & flag(:,1,1)~=0& resp_var(:,1,1,6)==2 ,1,1,9));
X3=abs(resp_var(ISI_fish'==5 & flag(:,1,1)~=0& resp_var(:,1,1,6)==3 ,1,1,9));
X10=abs(resp_var(ISI_fish'==5 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9));
Y1=abs(resp_var(ISI_fish'==5 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9));
%first side one and first side 2
kruskalwallis([X1;Y1],[ones(size(X1));2*ones(size(Y1))])
kruskalwallis([X2;Y1],[ones(size(X10));2*ones(size(Y1))])
kruskalwallis([X3;Y1],[ones(size(X10));2*ones(size(Y1))])

kruskalwallis([X1;X10],[ones(size(X1));2*ones(size(X10))])
kruskalwallis([X10;Y1],[ones(size(X10));2*ones(size(Y1))])

%compare ttc

%for 180
X1=(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,10));
X10=(resp_var(ISI_fish'==180 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,10));
Y1=(resp_var(ISI_fish'==180 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,10));
Y10=(resp_var(ISI_fish'==180 & flag(:,2,1)~=0& resp_var(:,2,1,6)==10,2,1,10));

%first side one and first side 2
kruskalwallis([X1;Y1],[ones(size(X1));2*ones(size(Y1))])
kruskalwallis([X1;X10],[ones(size(X1));2*ones(size(X10))])
kruskalwallis([X10;Y1],[ones(size(X10));2*ones(size(Y1))])
kruskalwallis([Y1;Y10],[ones(size(Y1));2*ones(size(Y10))])

%for 10

X1=(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,10));
X10=(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,10));
Y1=(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,10));
Y10=(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)==10,2,1,10));

%first side one and first side 2
kruskalwallis([X1;Y1],[ones(size(X1));2*ones(size(Y1))])
kruskalwallis([X1;X10],[ones(size(X1));2*ones(size(X10))])
kruskalwallis([X10;Y1],[ones(size(X10));2*ones(size(Y1))])
kruskalwallis([Y1;Y10],[ones(size(Y1));2*ones(size(Y10))])

% 
% 
% 
amps_10_3=[abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)<4 ,1,1,9));abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)>7,1,1,9));abs(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)<4,2,1,9))];
Group_10_3=[ones(size(abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)<4 ,1,1,9))));2*ones(size(abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)>7,1,1,9))));3*ones(size(abs(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)<4,2,1,9))))];
[P,ANOVATAB,STATS]=kruskalwallis(amps_10_3, Group_10_3);
figure
multcompare(STATS)

% %compare first three and last one and first on the other side

amps_10_3=[abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9));abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==2 ,1,1,9));abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==3 ,1,1,9));abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9));abs(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9))];
Group_10_3=[ones(size(abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==1 ,1,1,9))));2*ones(size(abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==2 ,1,1,9))));3*ones(size(abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==3 ,1,1,9))));4*ones(size(abs(resp_var(ISI_fish'==10 & flag(:,1,1)~=0& resp_var(:,1,1,6)==10,1,1,9))));5*ones(size(abs(resp_var(ISI_fish'==10 & flag(:,2,1)~=0& resp_var(:,2,1,6)==1,2,1,9))))];;
[P,ANOVATAB,STATS]=kruskalwallis(amps_10_3, Group_10_3);
figure
multcompare(STATS)


% 
% 
amps_5_3=[abs(resp_var(ISI_fish'==5 & flag(:,1,1)~=0& resp_var(:,1,1,6)<6 ,1,1,9));abs(resp_var(ISI_fish'==5 & flag(:,1,1)~=0& resp_var(:,1,1,6)>5,1,1,9));abs(resp_var(ISI_fish'==5 & flag(:,2,1)~=0& resp_var(:,2,1,6)<6,2,1,9))];
Group_5_3=[ones(size(abs(resp_var(ISI_fish'==5 & flag(:,1,1)~=0& resp_var(:,1,1,6)<6 ,1,1,9))));2*ones(size(abs(resp_var(ISI_fish'==5 & flag(:,1,1)~=0& resp_var(:,1,1,6)>5,1,1,9))));3*ones(size(abs(resp_var(ISI_fish'==5 & flag(:,2,1)~=0& resp_var(:,2,1,6)<6,2,1,9))))];
[P,ANOVATAB,STATS]=kruskalwallis(amps_5_3, Group_5_3);
figure
multcompare(STATS)
