%This is a function that calculates cluster types for neuronal data.
%Created by Haleh Fotowat
function [CCStimMot,trial_ave_mat,stim_mat] = ClusterTypes(repse,stim_type_cnt,F_dff_pool_all,ISI,Motor,fn,fps,per2p,stim_size_degrees)
%This is a function that clusters data through first dividing into cells
%with different tuning properties and then sorting individual cell types
%further based on their habituation dyanmics
Length_trial_type=(size(F_dff_pool_all,2)/stim_type_cnt);
%first calculate trail averages
trial_ave_mat=[];
stim_all=[];
%stim_mat=stim_size_degrees(1)*ones(stim_type_cnt,(fps+1)*repse);
f0=round(ISI/(2*per2p));

% t1=2;
% t2=t1+2;
t1=1;
t2=repse-1;
lb1=f0;
rb1=f0+fps;
for i=1:stim_type_cnt
    trial_ave_temp=0;
    for j=t1:t2
        lb=(j-1)*fps+f0+(i-1)*Length_trial_type;
        rb=j*fps+f0+(i-1)*Length_trial_type;
        trial_ave_temp=trial_ave_temp+(detrend((F_dff_pool_all(:,lb:rb))')');
    end
    trial_ave_mat=[trial_ave_mat,trial_ave_temp/(t2-t1+1)];
    stim_mat(i,(i-1)*(fps+1)+1:i*(fps+1))=(stim_size_degrees(lb1:rb1));
    stim_all=[stim_all,(stim_size_degrees(lb1:rb1))];
    %calculate peak amplitude for each trial

end
%for each cell, calculate the correlation coefficient with stimulus
%regressors for each stimulus type and with motor output (stim_type_cnt+1)
%values
CCStimMot=zeros(size(trial_ave_mat,1),stim_type_cnt+1);
for i=1:size(trial_ave_mat,1)
    for j=1:stim_type_cnt
        r=max(xcorr(trial_ave_mat(i,:),(stim_mat(j,:)-stim_mat(j,1)),3,'coeff'));
        CCStimMot(i,j)=r;
    end
    M=Motor{fn(i)}(1:size(F_dff_pool_all,2));
    Ml=M;
    Ml(M>0)=0;
    Mr=M;
    Mr(M<=0)=0;
    ra=max(xcorr(F_dff_pool_all(i,:),abs(M),20,'coeff'));%correlation with both left and right flick
    rl=max(xcorr(F_dff_pool_all(i,:),abs(Ml),20,'coeff'));%correlation with both left and right flick
    rr=max(xcorr(F_dff_pool_all(i,:),abs(Mr),20,'coeff'));
    CCStimMot(i,stim_type_cnt+1)=max([ra,rl,rr]);
end

    
end

