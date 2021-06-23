
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
    plot(state.time_myr,state.DEGASS,'r')
    plot(state.time_myr,state.BAS_AREA,'k')
    plot(state.time_myr,state.EVO,'g')
    plot(state.time_myr,state.W,'b')
    plot(state.time_myr,state.Bforcing,'m')
    %%%% Legend
    text(-590,2.4,'D','color','r')
    text(-590,2.2,'E','color','g')
    text(-590,2,'W','color','b')
    text(-590,1.8,'B','color','m')
    text(-590,1.6,'BAS','color','k')
    %%%% Title
    title('Forcings')

    %%% Corg fluxes
    subplot(4,5,2)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Flux (mol/yr)')
    %%%% plot this model
    plot(state.time_myr,state.mocb,'b')
    plot(state.time_myr,state.locb,'g')
    plot(state.time_myr,state.oxidw,'r')
    plot(state.time_myr,state.ocdeg,'k') 
    %%%% Legend
    text(-590,7e12,'mocb','color','b')
    text(-590,6e12,'locb','color','g')
    text(-590,5e12,'oxidw','color','r')
    text(-590,4e12,'ocdeg','color','k')
    %%%% Title
    title('C_{org} fluxes')

    %%% Ccarb fluxes
    subplot(4,5,3)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Flux (mol/yr)')
    %%%% plot this model
    plot(state.time_myr,state.silw,'r')
    plot(state.time_myr,state.carbw,'c')
    plot(state.time_myr,state.sfw,'b')
    plot(state.time_myr,state.mccb,'k') 
    %%%% Legend
    text(-590,33e12,'silw','color','r')
    text(-590,30e12,'carbw','color','c')
    text(-590,27e12,'sfw','color','b')
    text(-590,24e12,'mccb','color','k')
    %%%% Title
    title('C_{carb} fluxes')

    %%% S fluxes
    subplot(4,5,4)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    % ylim([0 5e12])
    ylabel('Fluxes (mol/yr)')
    %%%% plot this model
    plot(state.time_myr,state.mpsb,'k')
    plot(state.time_myr,state.mgsb,'c')
    plot(state.time_myr,state.pyrw,'r')
    plot(state.time_myr,state.pyrdeg,'m') 
    plot(state.time_myr,state.gypw,'b')
    plot(state.time_myr,state.gypdeg,'g') 
    %%%% Legend
    text(-590,1,9e12,'mpsb','color','k')
    text(-590,1.7e12,'mgsb','color','c')
    text(-590,1.5e12,'pyrw','color','r')
    text(-590,1.2e12,'pyrdeg','color','m')
    text(-590,1e12,'gypw','color','b')
    text(-590,0.8e12,'gypdeg','color','g')
    %%%% Title
    title('S fluxes')

    %%%% C SPECIES
    subplot(4,5,5)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Relative size')
    %%%% plot this model
    plot(state.time_myr,state.G/pars.G0,'k')
    plot(state.time_myr,state.C/pars.C0,'c')
    plot(state.time_myr,state.VEG,'g--')
    %%%% Legend
    text(-590,1.5,'VEG','color','g')
    text(-590,1.25,'G','color','k')
    text(-590,1,'C','color','c')
    %%%% Title
    title('C reservoirs')

    %%%% S SPECIES
    subplot(4,5,6)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('Relative size')
    %%%% plot this model
    plot(state.time_myr,state.PYR/pars.PYR0,'k')
    plot(state.time_myr,state.GYP/pars.GYP0,'c')
    %%%% Legend
    text(-590,1,'PYR','color','k')
    text(-590,0.9,'GYP','color','c')
    %%%% Title
    title('S reservoirs')

    %%% NUTRIENTS P N
    subplot(4,5,7)
    hold on
    box on
    xlim(pars.plotrange)
    ylim([0 3])
    xlabel('Time (Ma)')
    ylabel('Relative size')
    %%%% plot this model
    plot(state.time_myr,state.P/pars.P0,'b')
    plot(state.time_myr,state.N/pars.N0,'g')
    %%%% Legend
    text(-590,1.5,'P','color','b')
    text(-590,1,'N','color','g')
    %%%% Title
    title('Nutrient reservoirs')

    %%%% Forg and Fpy ratos
    subplot(4,5,8)
    hold on
    box on
    xlim(pars.plotrange)
    xlabel('Time (Ma)')
    ylabel('f_{org}, f_{py}')
    %%%% plot this model
    plot(state.time_myr,state.mocb ./ (state.mocb + state.mccb),'k')
    %%%% plot fpy
    plot(state.time_myr, state.mpsb ./ (state.mpsb + state.mgsb),'m')

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
    plot(state.time_myr,state.delta_G,'k')
    plot(state.time_myr,state.delta_C,'c')
    %%%% D34S other reservoirs (same fig)
    plot(state.time_myr,state.delta_PYR,'r--')
    plot(state.time_myr,state.delta_GYP,'g--')
    
    
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
