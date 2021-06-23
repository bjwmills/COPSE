%%%%%% COPSE for MATLAB
%%%%%%% ported by B Mills, 2013
%%%%%%% 2017 development version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run = COPSE_frontend(S)

    %%%%%%% remove structures from pervious runs 
    clear stepnumber
    clear pars
    clear forcings
    clear workingstate
    clear switches
    clear state
    clear rawoutput
    clear options
    clear geoldata
    clear rawoutput
    clear resample
    %%%%%%% set up global structures
    global stepnumber
    global pars
    global forcings
    global workingstate
    global sensanal
    global plotrun
    global sensparams
    %%%% global tuning variables
    global Gtune
    global Ctune
    global PYRtune
    global GYPtune
    global Atune
    global Otune
    global Stune
    
    %%%%%% check for sensitivity analysis
    if S >= 1
        sensanal = 1 ;
        plotrun = 0 ;
        pars.telltime = 0 ;
    else
        sensanal = 0 ;
        plotrun = 1 ;
        pars.telltime = 1 ;
    end
    
    %%%%%%% starting to load params
    if sensanal == 0 
        fprintf('setting parameters... \t')
        tic
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Flux values at present   %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% org C cycle
    pars.k_locb = 2.5e12 ;
    pars.k_mocb = 2.5e12 ;
    pars.k_ocdeg = 1.25e12 ;
    
    %%%% carb C cycle
    pars.k_ccdeg = 15e12 ;
    pars.k_carbw = 8e12 ;
    pars.k_sfw = 1.75e12 ;
    pars.k_mccb = pars.k_carbw + pars.k_ccdeg - pars.k_sfw ;
    pars.k_silw = pars.k_mccb - pars.k_carbw ;
    basfrac = 0.3 ;
    pars.k_granw = pars.k_silw * (1-basfrac) ;
    pars.k_basw = pars.k_silw * basfrac ;

    %%%% S cycle
    pars.k_mpsb = 0.7e12 ;
    pars.k_mgsb = 1.5e12 ;
    pars.k_pyrw = 4.5e11 ;
    pars.k_gypw = 1e12 ;
    pars.k_pyrdeg = 0.25e12 ; 
    pars.k_gypdeg = 0.5e12 ;
    %%%% P cycle
    pars.k_capb = 2e10 ;
    pars.k_fepb = 1e10 ;
    pars.k_mopb = 1e10 ;
    pars.k_phosw = 4.25e10 ;
    pars.k_landfrac = 0.0588 ;
    %%%% N cycle
    pars.k_nfix = 8.67e12 ;
    pars.k_denit = 4.3e12 ;
    
    %%%% reductant input
    pars.k_reductant_input = 0.4e12 ;  %%%% schopf and klein 1992
% %     pars.k_reductant_input = 0 ;

    %%%% fluxes calculated for steady state
    pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg - pars.k_reductant_input;

    %%%% Sr cycle
    pars.k_Sr_sedw = 17e9 ;
    pars.k_Sr_mantle = 7.3e9 ;
    pars.k_Sr_silw = 13e9 ;
    pars.k_Sr_granw = pars.k_Sr_silw * (1 - basfrac) ;
    pars.k_Sr_basw = pars.k_Sr_silw * basfrac ;
    pars.total_Sr_removal = pars.k_Sr_granw + pars.k_Sr_basw + pars.k_Sr_sedw + pars.k_Sr_mantle ;
    pars.k_Sr_sfw = pars.total_Sr_removal * ( pars.k_sfw / (pars.k_sfw + pars.k_mccb) ) ;
    pars.k_Sr_sedb = pars.total_Sr_removal * ( pars.k_mccb / (pars.k_sfw + pars.k_mccb) ) ;
    pars.k_Sr_metam = 13e9 ;

    %%%% others
    pars.k_oxfrac = 0.9975 ;
    Pconc0 = 2.2 ;
    Nconc0 = 30.9 ;
    pars.newp0 = 117 * min(Nconc0/16,Pconc0) ;
    %COPSE constant for calculating pO2 from normalised O2
    pars.copsek16 = 3.762 ;
    % R_PDB
    pars.rp = 0.0112372 ;
    pars.dmocb = -30 ;
    % R_VCDT
    pars.rv = 0.04416376 ;
    pars.dpyrb = -35 ;
    pars.DOC_res_0 = 3e18 * 40 ;
    

    %reservoir present day sizes (mol)
    pars.P0 = 3.1*10^15 ;
    pars.O0 = 3.7*10^19 ;
    pars.A0 = 3.193*10^18 ;
    pars.G0 = 1.25*10^21 ;
    pars.C0 = 5*10^21 ;
    pars.PYR0 = 1.8*10^20 ;
    pars.GYP0 = 2*10^20 ;
    pars.S0 = 4*10^19 ;
    pars.CAL0 = 1.397e19 ;
    pars.N0 = 4.35e16 ;
    pars.OSr0 = 1.2e17 ; %%% francois and walker 1992
    pars.SSr0 = 5e18 ;
    pars.U0 = 1.85e13 ; 

    %%%%%%% basic dependencies
    pars.a = 0.5 ; %oxidative weathering dependency on O2 concentration
    pars.b = 2 ; %marine organic carbon burial dependency on new production

    %%%%%%% weathering enhancement factor prior to vascular plant colonisation
    pars.plantenhance = 0.15 ;

    %%%%%%% transport limitation of weathering. 0 = fixed limit, 1 = rate
    %%%%%%% scales with global erosion/uplift rate
    pars.tlimitchoice = 0;

    %%%%%%% define maximum weathering rate relative to prsent day (k_transport)
    pars.TLIMIT = 10 ;

    %%fire feedback
    pars.kfire= 3 ;

    %%%% finished loading params
    if sensanal == 0 
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Forcings   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%% starting to load forcings
    if sensanal == 0 
        fprintf('loading forcings... \t')
        tic
    end

    %%%% load COPSE reloaded forcing set
    load( 'reloaded_forcings/forcings.mat' ) 

    %%%% load new forcings %%%%

    %%%% simple smoothing
    forcings.usmooth_2018 = xlsread('new_forcings/usmooth_2018.xlsx','','','basic') ;
    forcings.usmooth_2018(:,1) = forcings.usmooth_2018(:,1)*1e6 ; %%% correct Myr
    %%%% Degassing from Paleomap+VDM sbz length
    forcings.D_PALEOMAP_VDM = xlsread('new_forcings/D_paleomap_vdm_combined.xlsx','','','basic') ;
    forcings.D_PALEOMAP_VDM(:,1) = forcings.D_PALEOMAP_VDM(:,1)*1e6 ; %%% correct Myr
    %%%% extended granite area
    forcings.GA_revised = xlsread('new_forcings/GA_revised.xlsx','','','basic') ;
    forcings.GA_revised(:,1) = forcings.GA_revised(:,1)*1e6 ; %%% correct Myr
    %%%% new degassing
    forcings.D_SBZ_RIFT = xlsread('new_forcings/D_sbz_rift.xlsx','','','basic') ;
    forcings.D_SBZ_RIFT(:,1) = forcings.D_SBZ_RIFT(:,1)*1e6 ; %%% correct Myr
    %%%% PG with cryo increase
    forcings.cryo_PG = xlsread('new_forcings/cryo_PG.xlsx','','','basic') ;
    forcings.cryo_PG(:,1) = forcings.cryo_PG(:,1)*1e6 ; %%% correct Myr
    %%%% new BA including the altered degassing forcing
    forcings.GR_BA = xlsread('new_forcings/GR_BA.xlsx','','','basic') ;
    forcings.GR_BA(:,1) = forcings.GR_BA(:,1)*1e6 ; %%% correct Myr

    %%%%% finished loading forcings
    if sensanal == 0 
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Sensitivity analysis   %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if sensanal == 1
        %%%% generate random number in [-1 +1]
        randminusplus1 = 2*(0.5 - rand) ;
        randminusplus2 = 2*(0.5 - rand) ;
        randminusplus3 = 2*(0.5 - rand) ;
        randminusplus4 = 2*(0.5 - rand) ;

        %%%%%%% parameter space to test
        sensparams.DEGASS = 1 + 0.2*randminusplus1 ;
        sensparams.UPLIFT = 1 + 0.2*randminusplus2 ;
        sensparams.BAS_AREA = 1 + 0.2*randminusplus3 ;
        sensparams.GRAN_AREA = 1 + 0.2*randminusplus4 ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Initialise   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% run beginning
    if sensanal == 0 
        fprintf('Beginning run: \n')
    end
    %%%%%%% model timeframe in years (0 = present day)
    pars.whenstart = - 900e6 ;
%     pars.whenstart = - 600e6 ;
    pars.whenend = 0 ;

    % pars.plotrange = [-620 -590] ;
    pars.plotrange = [-800 0] ;

    %%%%%%% set number of model steps to take before beiling out
    pars.bailnumber = 1e5;

    %%%%%%% display every n model steps whilst running 
    pars.display_resolution = 200 ;

    %%%%%%% set maximum step size for solver
    options = odeset('maxstep',1e6) ;

    %%%% set stepnumber to 1
    stepnumber = 1 ;

    %%%%%%% set starting reservoir sizes 
    pars.pstart = pars.P0*2.5;
    pars.tempstart = 288;
    pars.CAL_start = pars.CAL0;
    pars.N_start = pars.N0;
    pars.OSr_start = pars.OSr0;
    pars.SSr_start = pars.SSr0;
    pars.U_start = pars.U0 ;
    pars.delta_A_start = 0 ;
    pars.delta_S_start = 20 ;
    pars.delta_G_start = -27 ;
    pars.delta_C_start = -2 ;
    pars.delta_PYR_start = -5 ;
    pars.delta_GYP_start = 20 ;
    pars.delta_OSr_start = 0.708 ;
    pars.delta_SSr_start = 0.708 ;
    pars.delta_U_start = -0.4 ;
    pars.DOC_start = pars.DOC_res_0 ;
%     pars.DOC_start = 0 ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Initial parameter tuning  %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(Gtune) == 0
        pars.ostart = pars.O0 * abs( Otune )  ;
        pars.astart = pars.A0 * abs( Atune ) ;
        pars.sstart = pars.S0 * abs( Stune ) ;
        pars.gstart = pars.G0 * abs( Gtune ) ;
        pars.cstart = pars.C0 * abs( Ctune ) ;
        pars.pyrstart = pars.PYR0 * abs( PYRtune ) ;
        pars.gypstart = pars.GYP0 * abs( GYPtune ) ; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% if no tuning use previously tuned values
    if isempty(Gtune) == 1

%     outputs = [ 0.5 1.2 2.5 1 0.1 1 3] ;
    outputs = [ 0.33 1 1.5 0.5 0.1 1 3] ;

        pars.gstart = pars.G0 * outputs(1) ;
        pars.cstart = pars.C0 * outputs(2) ;
        pars.pyrstart = pars.PYR0 * outputs(3) ;
        pars.gypstart = pars.GYP0 * outputs(4) ; 
        pars.ostart = pars.O0 * outputs(5)  ;
        pars.sstart = pars.S0 * outputs(6) ;
        pars.astart = pars.A0 * outputs(7) ;

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%% model start state
    pars.startstate(1) = pars.pstart ;
    pars.startstate(2) = pars.ostart ;
    pars.startstate(3) = pars.astart ;
    pars.startstate(4) = pars.sstart ;
    pars.startstate(5) = pars.gstart ;
    pars.startstate(6) = pars.cstart ;
    pars.startstate(7) = pars.pyrstart ;
    pars.startstate(8) = pars.gypstart ;
    pars.startstate(9) = pars.tempstart ;
    pars.startstate(10) = pars.CAL_start ;
    pars.startstate(11) = pars.N_start ;
    pars.startstate(12) = pars.gstart * pars.delta_G_start ;
    pars.startstate(13) = pars.cstart * pars.delta_C_start ;
    pars.startstate(14) = pars.pyrstart * pars.delta_PYR_start ;
    pars.startstate(15) = pars.gypstart * pars.delta_GYP_start ;
    pars.startstate(16) = pars.astart * pars.delta_A_start ;
    pars.startstate(17) = pars.sstart * pars.delta_S_start ;
    pars.startstate(18) = pars.OSr_start ;
    pars.startstate(19) = pars.OSr_start * pars.delta_OSr_start ;
    pars.startstate(20) = pars.SSr_start ;
    pars.startstate(21) = pars.SSr_start * pars.delta_SSr_start ;
    pars.startstate(22) = pars.DOC_start ;
    pars.startstate(23) = pars.U_start ;
    pars.startstate(24) = pars.U_start * pars.delta_U_start ;

    %%%% note model start time
    tic

    %%%%%%% run the system 
    [rawoutput.T,rawoutput.Y] = ode15s(@COPSE_equations,[pars.whenstart pars.whenend],pars.startstate,options);





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    %%%% size of output 
    pars.output_length = length(rawoutput.T) ;


    if sensanal == 0

        %%%%%%%%%% model finished output to screen
        fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length ) 
        toc
        
        %%%%%% check conservation of C, S 
        Cano = 100 * ( workingstate.res_C(end) - workingstate.res_C(1) )  ;
        Sano = 100 * ( workingstate.res_S(end) - workingstate.res_S(1) )  ;       
        %%%%%% check conservation of C, S isotopes
        Cisoano = 100 * ( workingstate.iso_res_C(end) - workingstate.iso_res_C(1) ) / workingstate.iso_res_C(1) ;
        Sisoano = 100 * ( workingstate.iso_res_S(end) - workingstate.iso_res_S(1) ) / workingstate.iso_res_S(1) ;       
        fprintf('C anomoly (mol): %e \n', Cisoano  ) 
        fprintf('S anomoly (mol): %e \n', Sisoano ) 
        fprintf('C*isotope anomoly (percent): %e \n', Cisoano  ) 
        fprintf('S*isotope anomoly (percent): %e \n', Sisoano ) 

    end



    %%%%%%%%% print final model states using final state for each timepoint
    %%%%%%%%% during integration
    
    if sensanal == 0
    fprintf('assembling state vectors... \t')
    tic
    end
    
    %%%% trecords is index of shared values between ode15s output T vector and
    %%%% model recorded workingstate t vector
    [sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

    %%%%%% assemble output state vectors
    field_names = fieldnames(workingstate) ;
    for numfields = 1:length(field_names)
        eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
    end

    %%%%%% save state
    run.state = state ;
    run.pars = pars ;
    run.forcings = forcings ;
    
    if sensanal == 0
        %%%%%% done message
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% only plot if no tuning structure exists
    if isempty(Gtune) == 1
        if plotrun == 1
            COPSE_plot
        end
    end

    

end