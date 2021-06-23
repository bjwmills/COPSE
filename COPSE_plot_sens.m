
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Plotting script   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% define colours
c_mean = [0.2 0.6 0.6] ;
c_std = [0.3 0.7 0.7] ;
c_range = [ 0.4 0.8 0.8] ;

%%%% output to screen
fprintf('running sens plotting script... \t')
tic

%%%%%%% make figure
figure('Color',[0.80 0.80 0.70])

%%%% load geochem data
load('data/data.mat')

%%%% make column vector
sens.time_myr = sens.time_myr(:,1) ;

%%%% d13C record
subplot(3,3,1)
hold on
box on
xlabel('Time (Ma)')
ylabel('\delta^{13}C_{carb}')
%%%% plot data comparison
plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
%%%% plot this model
plot((sens.time_myr),mean(sens.delta_mccb,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.delta_mccb,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.delta_mccb,[],2),'linewidth',0.5,'color',c_range)

%%%% d34S record
subplot(3,3,2)
hold on
box on
xlabel('Time (Ma)')
ylabel('\delta^{34}S_{sw}')
%%%% plot data comparison
plot(d34s_x,d34s_y,'.-','color',[0.8 0.8 0.8])
%%%% plot this model
plot((sens.time_myr),mean(sens.d34s_S,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.d34s_S,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.d34s_S,[],2),'linewidth',0.5,'color',c_range)

%%%% Ocean 87Sr/86Sr 
subplot(3,3,3)
hold on
box on
ylim([0.706 0.71])
xlabel('Time (Ma)')
ylabel('^{87}Sr/^{86}Sr seawater')
%%%% plot data comparison
plot(sr_x,sr_y,'.-','color',[0.8 0.8 0.8])
%%%% plot this model
plot((sens.time_myr),mean(sens.delta_OSr,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.delta_OSr,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.delta_OSr,[],2),'linewidth',0.5,'color',c_range)

%%%% SO4
subplot(3,3,4)
hold on
box on
xlabel('Time (Ma)')
ylabel('Marine SO_{4} (mM)')
%%%% plot this model
plot((sens.time_myr),mean(sens.SmM,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.SmM,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.SmM,[],2),'linewidth',0.5,'color',c_range)

%%%% O2 (%) 
subplot(3,3,5)
hold on
box on
xlabel('Time (Ma)')
ylabel('Atmospheric O_{2} (%)')
%%%% plot this model
plot((sens.time_myr),mean(sens.mrO2./0.21,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.mrO2./0.21,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.mrO2./0.21,[],2),'linewidth',0.5,'color',c_range)


%%%% CO2ppm
subplot(3,3,6)
hold on
box on
ylim([1.5 4.5])
yticks([1 2 3 4])
yticklabels({'10','100','1000','10,000'})
xlabel('Time (Ma)')
ylabel('Atmospheric CO_{2} (ppm)')
%%%% plot this model
plot((sens.time_myr),log10(mean(sens.CO2ppm,2)),'linewidth',1,'color',c_mean)
plot((sens.time_myr),log10(max(sens.CO2ppm,[],2)),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),log10(min(sens.CO2ppm,[],2)),'linewidth',0.5,'color',c_range)

%%%% TEMP
subplot(3,3,7)
hold on
box on
ylim([5 40])
xlabel('Time (Ma)')
ylabel('GAST (C)')
%%%% plot this model
plot((sens.time_myr),mean(sens.T_gast,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.T_gast,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.T_gast,[],2),'linewidth',0.5,'color',c_range)


%%%% DOC res
subplot(3,3,8)
hold on
box on
xlabel('Time (Ma)')
ylabel('DOC res (%)')
%%%% plot this model
plot((sens.time_myr),mean(sens.DOC_res,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.DOC_res,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.DOC_res,[],2),'linewidth',0.5,'color',c_range)


%%%% ANOX
subplot(3,3,9)
hold on
box on
xlabel('Time (Ma)')
ylabel('ANOX')
%%%% plot this model
plot((sens.time_myr),mean(sens.ANOX,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.ANOX,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.ANOX,[],2),'linewidth',0.5,'color',c_range)




%%%%% plotting script finished
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Cleanup workspace   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear stepnumber
clear u
clear numfields
clear trecords
clear finalrecord
clear field_names
clear n
clear veclength
clear xvec
clear yvec
clear endtime


%%%%% plotting script finished
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )
