
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Plotting script   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% only plot if no tuning structure exists
if isempty(Gtune) == 1

    %%%% output to screen
    fprintf('running plotting script... \t')
    tic


    %%%%%%% make figure
    figure('Color',[0.80 0.80 0.70])
    
    %%%% load geochem data
    load('data/data.mat')


    %%%% GLOBAL FORCINGS
    subplot(4,5,1)
    hold on
    box on
    xlim(pars.plotrange)
    ylim([0 2.5])
    xlabel('Time (Ma)')
    ylabel('Relative forcing')
    %%%% plot this model
    plot(state.time_myr,state.DEGASS,'r','displayname','D')
    plot(state.time_myr,state.BAS_AREA,'k','displayname','BA')
     plot(state.time_myr,state.GRAN_AREA,'k--','displayname','GA')
    plot(state.time_myr,state.EVO,'g','displayname','E')
    plot(state.time_myr,state.W,'b','displayname','W')
    plot(state.time_myr,state.Bforcing,'m','displayname','B')
    %%%% Title
    title('Forcings')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')

    %%% Corg fluxes
    subplot(4,5,2)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Flux (mol/yr)')
    %%%% plot this model
    plot(state.time_myr,state.mocb,'b','displayname','mocb')
    plot(state.time_myr,state.locb,'g','displayname','locb')
    plot(state.time_myr,state.oxidw,'r','displayname','oxidw')
    plot(state.time_myr,state.ocdeg,'k','displayname','ocdeg') 
    %%%% Title
    title('C_{org} fluxes')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')

    %%% Ccarb fluxes
    subplot(4,5,3)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Flux (mol/yr)')
    %%%% plot this model
    plot(state.time_myr,state.silw,'r','displayname','silw')
    plot(state.time_myr,state.carbw,'c','displayname','carbw')
    plot(state.time_myr,state.sfw,'b','displayname','sfw')
    plot(state.time_myr,state.mccb,'k','displayname','mccb') 
    %%%% Title
    title('C_{carb} fluxes')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')
   
    %%% S fluxes
    subplot(4,5,4)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    % ylim([0 5e12])
    ylabel('Fluxes (mol/yr)')
    %%%% plot this model
    plot(state.time_myr,state.mpsb,'k','displayname','mpsb')
    plot(state.time_myr,state.mgsb,'c','displayname','mgsb')
    plot(state.time_myr,state.pyrw,'r','displayname','pyrw')
    plot(state.time_myr,state.pyrdeg,'m','displayname','pyrdeg') 
    plot(state.time_myr,state.gypw,'b','displayname','gypw')
    plot(state.time_myr,state.gypdeg,'g','displayname','gypdeg') 
    %%%% Title
    title('S fluxes')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')
    
    %%%% C SPECIES
    subplot(4,5,5)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Relative size')
    %%%% plot this model
    plot(state.time_myr,state.G/pars.G0,'k','displayname','G')
    plot(state.time_myr,state.C/pars.C0,'c','displayname','C')
    plot(state.time_myr,state.VEG,'g--','displayname','VEG')
    %%%% Title
    title('C reservoirs')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')
    
    %%%% S SPECIES
    subplot(4,5,6)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Relative size')
    %%%% plot this model
    plot(state.time_myr,state.PYR/pars.PYR0,'k','displayname','PYR')
    plot(state.time_myr,state.GYP/pars.GYP0,'c','displayname','GYP')
    %%%% Title
    title('S reservoirs')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')
    
    %%% NUTRIENTS P N
    subplot(4,5,7)
    hold on
    box on
    xlim(pars.plotrange)
    ylim([0 3])
    xlabel('Time (Ma)')
    ylabel('Relative size')
    %%%% plot this model
    plot(state.time_myr,state.P/pars.P0,'b','displayname','P')
    plot(state.time_myr,state.N/pars.N0,'g','displayname','N')
    %%%% Title
    title('Nutrient reservoirs')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')
    
    %%%% Forg and Fpy ratos
    subplot(4,5,8)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('f_{org}, f_{py}')
    %%%% plot this model
    plot(state.time_myr,state.mocb ./ (state.mocb + state.mccb),'k','displayname','f_{org}')
    %%%% plot fpy
    plot(state.time_myr, state.mpsb ./ (state.mpsb + state.mgsb),'m','displayname','f_{py}')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')
    
    %%%% d13C record
    subplot(4,5,9)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('\delta^{13}C_{carb}')
    %%%% plot data comparison
    plot(d13c_x,d13c_y,'.-','color',[0.8 0.8 0.8])
    %%%% plot this model
    plot(state.time_myr,state.delta_mccb,'k')

    %%%% d34S record
    subplot(4,5,10)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('\delta^{34}S_{sw}')
    %%%% plot data comparison
    plot(d34s_x,d34s_y,'.-','color',[0.8 0.8 0.8])
    %%%% plot this model
    plot(state.time_myr,state.d34s_S,'k')

    %%%% Ocean 87Sr/86Sr 
    subplot(4,5,11)
    hold on
    box on
    xlim(pars.plotrange)
    ylim([0.706 0.71])
    xlabel('Time (Ma)')
    ylabel('^{87}Sr/^{86}Sr seawater')
    %%%% plot data comparison
    plot(sr_x,sr_y,'.-','color',[0.8 0.8 0.8])
    %%%% plot this model
    plot(state.time_myr,state.delta_OSr,'k')

    %%%% SO4
    subplot(4,5,12)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Marine SO_{4} (mM)')
    %%%% plot this model
    plot(state.time_myr,(state.S./pars.S0)*28,'k')

    %%%% O2 (%) 
    subplot(4,5,13)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Atmospheric O_{2} (%)')
    %%%% plot this model
    plot(state.time_myr,state.mrO2.*100,'k')

    %%%% CO2ppm
    subplot(4,5,14)
    hold on
    box on
    xlim(pars.plotrange)
    ylim([1.5 4.5])
    yticks([1 2 3 4])
    yticklabels({'10','100','1000','10,000'})
    xlabel('Time (Ma)')
    ylabel('Atmospheric CO_{2} (ppm)')
    %%%% plot this model
    plot(state.time_myr,log10(state.RCO2.*280),'k')

    %%%% TEMP
    subplot(4,5,15)
    hold on
    box on
    xlim(pars.plotrange)
    ylim([5 40])
    xlabel('Time (Ma)')
    ylabel('GAST (C)')
    %%%% plot this model
    plot(state.time_myr,state.tempC,'k')

    
    %%%% DOC res
    subplot(4,5,16)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('DOC res (%)')
    %%%% plot this model
    plot(state.time_myr,state.DOC_res,'k')
   
    %%%% DOC oxidation
    subplot(4,5,17)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('DOC ox flux (mol/yr)')
    %%%% plot this model
    plot(state.time_myr,state.DOC_ox,'k')
   

    %%%% ANOX
    subplot(4,5,18)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('ANOX')
    %%%% plot this model
    plot(state.time_myr,state.ANOX,'k')
    
    
%     %%%% 
%     subplot(4,5,19)
%     hold on
%     box on
%     xlim(pars.plotrange)
%     xlabel('Time (Ma)')
%     ylabel('[U]')
%     %%%% plot this model
%     plot(state.time_myr,state.U,'k')

    %%%% D13C other reservoirs
    subplot(4,5,19)
    hold on
    box on
    xlim([pars.whenstart/1e6 pars.whenend/1e6])
    xlabel('Time (Ma)')
    ylabel('d13c C and S crust')
    %%%% plot this model
    plot(state.time_myr,state.delta_G,'k','displayname','deltaG')
    plot(state.time_myr,state.delta_C,'c','displayname','deltaC')
    %%%% D34S other reservoirs (same fig)
    plot(state.time_myr,state.delta_PYR,'r--','displayname','deltaPYR')
    plot(state.time_myr,state.delta_GYP,'g--','displayname','deltaGYP')
    %%%% legend
    l = legend ;
    set(l,'fontsize',6)
    set(l,'edgecolor','none')
    set(l,'location','northwest')
    
    %%%% 
    subplot(4,5,20)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('\delta^{238}U')
    %%%% plot this model
    plot(state.time_myr,state.d238U_sw,'k')
    

    %%%%% plotting script finished
    fprintf('Done: ')
    endtime = toc ;
    fprintf('time (s): %d \n', endtime )


end





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
