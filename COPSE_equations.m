function dy = COPSE_equations(t,y)

%%%% COPSE V2.1 (Carbon Oxygen Phosphorus Sulfur Evolution)
%%%% As used in Tostevin and Mills (2020) Interface Focus
%%%% Coded by Benjamin JW Mills // b.mills@leeds.ac.uk


%%%%%%% setup dy array
dy = zeros(21,1);  

%%%%%%% set up global parameters
global stepnumber
global pars
global forcings
global workingstate
global sensanal
global sensparams


%%%%%%%%%%%%% get variables from Y to make working easier
P = y(1) ;
O = y(2) ;
A = y(3) ;
S = y(4) ;
G = y(5) ;
C = y(6) ;
PYR = y(7) ;
GYP = y(8) ;
% TEMP = y(9);
% CAL = y(10) ;
N = y(11) ;
OSr = y(18) ;
SSr = y(20) ;
dSSr = y(21)/y(20) ;
DOC_res = y(22) ;
U = y(23) ;
d238U_sw = y(24) / y(23) ;

%%%% geological time in Ma
t_geol = t*(1e-6) ;

%%%%%%% calculate isotopic fractionation of reservoirs
delta_G = y(12)/y(5);
delta_C = y(13)/y(6);
delta_GYP  = y(15)/y(8);
delta_PYR  = y(14)/y(7);

%%%%%%% atmospheric fraction of total CO2, atfrac(A)
atfrac0 = 0.01614 ;
%%%%%%% constant
% atfrac = 0.01614 ;
%%%%%%% variable
atfrac = atfrac0 * (A/pars.A0) ;

%%%%%%%% calculations for pCO2, pO2
RCO2 = (A/pars.A0)*(atfrac/atfrac0) ;
CO2atm = RCO2*(280e-6) ;

%%%%% mixing ratio of oxygen (not proportional to O reservoir)
mrO2 = ( O/pars.O0 )  /   ( (O/pars.O0)  + pars.copsek16 ) ;
%%%%% relative moles of oxygen 
RO2 =  O/pars.O0 ;
%%%%% pO2 = mixing ratio * atmospheric pressure



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Interpolate forcings for this timestep   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% COPSE Reloaded forcing set
CP_reloaded = interp1qr( 1e6 * forcings.t', forcings.CP' , t ) ;
E_reloaded = interp1qr( 1e6 * forcings.t', forcings.E' , t ) ;
W_reloaded = interp1qr( 1e6 * forcings.t', forcings.W' , t ) ;
coal_reloaded = interp1qr( 1e6 * forcings.t', forcings.coal' , t ) ;

%%%% Additional forcings
GA_revised = interp1qr( forcings.GA_revised(:,1) , forcings.GA_revised(:,2) , t ) ;
D_sbz_rift = interp1qr( forcings.D_SBZ_RIFT(:,1) , forcings.D_SBZ_RIFT(:,2) , t ) ;
U_smooth2018 = interp1qr( forcings.usmooth_2018(:,1) , forcings.usmooth_2018(:,2)./forcings.usmooth_2018(end,2) , t ) ;
cryo_PG = interp1qr( forcings.cryo_PG(:,1) , forcings.cryo_PG(:,2) , t ) ;
GR_BA = interp1qr( forcings.GR_BA(:,1) , forcings.GR_BA(:,2) , t ) ;
epsilon_new = interp1qr([-1000e6 -466e6 -444e6 -411e6 -399e6 0 ]',[0.75 0.75 1.5 1.5 1 1 ]',t) ;

%%%% S inputs
GYP_INPUT = interp1qr([-1000e6 -581e6 -580e6 -571e6 -570e6 0]',[0 0 7 7 0 0]',t) ;
PYR_INPUT = interp1qr([-1000e6 -581e6 -580e6 -571e6 -570e6 0]',[0 0 7 7 0 0]',t) ;
pyrburialfrac = 0.8 ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Choose forcing functions  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UPLIFT = U_smooth2018 ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEGASS = D_reloaded ;
DEGASS = D_sbz_rift ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W = 1 ;
W = W_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVO = 1 ;
EVO = E_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPland_relative = 1 ;
CPland_relative = CP_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPSILON = 1 ;
EPSILON = epsilon_new ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bforcing = interp1qr([-1000 -150 -100 0]',[0.75 0.75 1 1]',t_geol) ;
% Bforcing = B_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PG = 1 ;
PG = cryo_PG ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAS_AREA = 1 ;
BAS_AREA = GR_BA ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CARB_AREA = 1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAN_AREA = 1 ;
GRAN_AREA = GA_revised ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ORG_AREA = 1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COALF = 1 ;
COALF = coal_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CAL = 1 ;
% CAL = Ca_reloaded ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEOG = GEOG_royer ;
GEOG = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% weathering relationships
silconst = 0.33 ;
carbconst = 0.9 ;
% silconst = 1 ;
% carbconst = 1 ;

%%%% reductant input
REDUCT_FORCE = DEGASS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Sensitivity analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sensanal == 1
    DEGASS = DEGASS * sensparams.DEGASS ;
    UPLIFT = UPLIFT * sensparams.UPLIFT ;
    BAS_AREA = BAS_AREA * sensparams.BAS_AREA ;
    GRAN_AREA = GRAN_AREA * sensparams.GRAN_AREA ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Calculate variables   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Berner temperature for COPSE reloaded, climate sensitivity variable
%%%% use Earth system sensitivity 
climsens = 5 ; 
%%%% temperature for weathering takes into account land surface temperature
%%%% GEOG
TEMP_land = 288 + climsens*(log(RCO2)/log(2)) - 7.4*(t_geol/-570) + GEOG ; 
%%%% Global average surface tmeperature doesn't take into acocunt GEOG
TEMP_gast = 288 + climsens*(log(RCO2)/log(2)) - 7.4*(t_geol/-570) ; 

%%%% relcalculate low lat temp for surface processes 
tgrad = 0.66 ;
% tgrad = 1 ;
tc = 288*(1-tgrad) ;
Tsurf = TEMP_land*tgrad + tc + 10 ; %%% low lat temp 25C at present
% Tsurf = TEMP ;

%%%% effect of temp on VEG 
V_T = 1 - (( (Tsurf - 298)/25 )^2) ;

%%%% effect of CO2 on VEG
P_atm = CO2atm*1e6 ;
P_half = 183.6 ;
P_min = 10 ;
V_co2 = (P_atm - P_min) / (P_half + P_atm - P_min) ;

%%%% effect of O2 on VEG
V_o2 = 1.5 - 0.5*(O/pars.O0) ; %% COPSE forumla

%%%% full VEG limitation
V_npp = 2*EVO*V_T*V_o2*V_co2 ;

%%% fire feedback
% ignit = max(586.2*(mrO2)-122.102 , 0  ) ;
%%%% COPSE reloaded
ignit = min(max(48*mrO2 - 9.08 , 0) , 5 ) ;
firef = pars.kfire/(pars.kfire - 1 + ignit) ;

%%% Mass of terrestrial biosphere
VEG = V_npp * firef ;

%%%%%% basalt and granite temp dependency - direct and runoff
f_T_bas =  exp(0.0608*(Tsurf-298)) * ( (1 + 0.038*(Tsurf - 298))^0.65 ) ; %%% 42KJ/mol
f_T_gran =  exp(0.0724*(Tsurf-298)) * ( (1 + 0.038*(Tsurf - 298))^0.65 ) ; %%% 50 KJ/mol
g_T = 1 + 0.087*(Tsurf - 298) ;

%%%% COPSE reloaded fbiota
V = VEG ;
f_biota = ( 1 - min( V*W , 1 ) ) * pars.plantenhance * (RCO2^0.5) + (V*W) ;

%%%% basalt and granite weathering
basw = pars.k_basw * BAS_AREA * PG * f_biota * f_T_bas ;
granw = pars.k_granw * UPLIFT^silconst * GRAN_AREA * PG * f_biota * f_T_gran ;
%%% silicate weathering
silw = basw + granw ;

%%%% carbonate weathering 
carbw = pars.k_carbw * CARB_AREA * UPLIFT^carbconst * PG * f_biota * g_T * (C/pars.C0) ;
% carbw = pars.k_carbw * CARB_AREA * UPLIFT * PG * f_biota * g_T ;

%%%% oxidative weathering 
oxidw = pars.k_oxidw*UPLIFT^silconst*ORG_AREA*(G/pars.G0)*((O/pars.O0)^pars.a) ;

%%% pyrite weathering
pyrw = pars.k_pyrw*UPLIFT^silconst*(PYR/pars.PYR0)  ;

%%% gypsum weathering 
gypw = pars.k_gypw*(GYP/pars.GYP0)*(carbw/pars.k_carbw) ;

%%% evap dissolution and pyrite delivery
evapdis = pars.k_gypw*GYP_INPUT ;
pulse_pyr = pars.k_pyrw*PYR_INPUT ;

%%%%% seafloor weathering, revised following Brady and Gislason but not directly linking to CO2
f_T_sfw = exp(0.0608*(TEMP_gast-288)) ; 
sfw = pars.k_sfw * f_T_sfw * DEGASS ; %%% assume spreading rate follows degassing here

%%%% carbonate burial via alkalinity balance
mccb = carbw + silw ;

%%%%%%% Degassing 
ocdeg = pars.k_ocdeg*DEGASS*(G/pars.G0) ;
ccdeg = pars.k_ccdeg*DEGASS*(C/pars.C0)*Bforcing ;
pyrdeg = pars.k_pyrdeg*(PYR/pars.PYR0)*DEGASS;
gypdeg = pars.k_gypdeg*(GYP/pars.GYP0)*DEGASS;

%%%% COPSE reloaded P weathering
pfrac_silw = 0.8 ;
pfrac_carbw = 0.14 ;
pfrac_oxidw = 0.06 ;
%%%% p weathering modified for young and ancient reservoirs
%%%% COPSE reloaded formula
phosw = EPSILON * pars.k_phosw*( (pfrac_silw)*( silw/pars.k_silw )  +   (pfrac_carbw)*( carbw/pars.k_carbw ) +  (pfrac_oxidw)*(  oxidw/ pars.k_oxidw )  )  ;

%%%% COPSE reloaded
k_aq = 0.8 ;
pland = pars.k_landfrac * VEG * phosw * ( k_aq + ( 1 - k_aq )*COALF ) ;
pland0 = pars.k_landfrac*pars.k_phosw;
psea = phosw - pland ;

%%%% convert total reservoir moles to micromoles/kg concentration    
Pconc = ( P/pars.P0 ) * 2.2 ;
Nconc = ( N/pars.N0 ) * 30.9 ;
newp = 117 * min(Nconc/16,Pconc) ;    

%%%%%% OCEAN ANOXIC FRACTION
k_anox = 10 ; 
k_u = 0.4 ;
ANOX = 1 / ( 1 + exp( -1 * k_anox * ( k_u * (newp/pars.newp0) - (O/pars.O0) ) ) ) ;

%%%% bioturbation forcing
f_biot = interp1qr([-1000e6 -525e6 -520e6 0]',[0 0 1 1]',t);
CB = interp1qr([0 1]',[1.2 1]',f_biot) ;

%%%%% carbon burial
mocb = pars.k_mocb*((newp/pars.newp0)^pars.b) * CB ;
locb = pars.k_locb*(pland/pland0)*CPland_relative  ;

% PYR burial function (COPSE)
fox= 1/(O/pars.O0) ; 

%%%% sulfur burial
mpsb = pars.k_mpsb*(S/pars.S0)*fox*(mocb/pars.k_mocb)+ pyrburialfrac*(evapdis + pulse_pyr) ;
mgsb = pars.k_mgsb*(S/pars.S0) + (1-pyrburialfrac)*(evapdis + pulse_pyr) ;

CNsea = 37.5 ;
% mopb = ( mocb/CPsea ) ;
monb = mocb/CNsea ;

%%%% P burial with bioturbation on
CPbiot = 250 ;
CPlam = 1000 ;
mopb = mocb*( (f_biot/CPbiot) + ( (1-f_biot)/CPlam ) ) ;

capb = pars.k_capb*( mocb/pars.k_mocb ) ;
% capb = pars.k_capb*( (newp/pars.newp0)^pars.b );

%%%% reloaded
fepb = (pars.k_fepb/pars.k_oxfrac)*(1-ANOX)*(P/pars.P0) ;

%%%%% nitrogen cycle
%%%% COPSE reloaded
if (N/16) < P
    nfix = pars.k_nfix * ( ( ( P - (N/16)  ) / (  pars.P0 - (pars.N0/16)    ) )^2 ) ;
else
    nfix = 0 ;
end

denit = pars.k_denit * ( 1 + ( ANOX / (1-pars.k_oxfrac) )  ) * (N/pars.N0) ;




%%%% DOC oxidation
kbig = 1e14 ;
if DOC_res > 1e12
    DOC_ox = sigmf((1-ANOX),[300 0.5]) .* kbig *(DOC_res/pars.DOC_res_0) ;   
else
    DOC_ox = 0 ;
end

%%%% reductant input
reductant_input = pars.k_reductant_input * REDUCT_FORCE ;
% reductant_input = pars.k_reductant_input * DEGASS ;

Gsub = 0 ;
Csub = 0 ;
PYRsub = 0 ;
GYPsub = 0 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Reservoir calculations  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Phosphate
CPDOC = 250 ;
dy(1) = psea - mopb - capb - fepb + ( DOC_ox / CPDOC ) ;
% dy(1) = psea - mopb - capb - fepb  ;

%%% Oxygen
dy(2) = locb + mocb - oxidw  - ocdeg   + 2*(mpsb - pyrw  - pyrdeg - pulse_pyr) - DOC_ox - reductant_input ;

%%% Carbon dioxide
dy(3) = -locb - mocb + oxidw + ocdeg + ccdeg + carbw - mccb - sfw + DOC_ox + reductant_input ;

%%% Sulphate
dy(4) = gypw + pyrw - mgsb - mpsb + gypdeg + pyrdeg  + evapdis + pulse_pyr ;

%%%Buried organic C
dy(5) = locb + mocb - oxidw - ocdeg - Gsub;
% dy(5) = 0 ;

%%% Buried carb C 
dy(6) = mccb + sfw - carbw - ccdeg - Csub;

%%% Buried pyrite S
dy(7) = mpsb - pyrw - pyrdeg - pulse_pyr - PYRsub;

%%% Buried gypsum S 
dy(8) = mgsb - gypw - gypdeg - evapdis - GYPsub ;

%%%% Nitrate
dy(11) = nfix - denit - monb;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Isotope reservoirs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% d13c and d34s for forwards model
d13c_A = y(16) / y(3) ;
d34s_S = y(17) / y(4) ;

%%%% carbonate fractionation
delta_o = atfrac*( (9483/TEMP_gast)  - 23.89 ) ; 
d_mccb = delta_o + 15.1 - (4232/TEMP_gast) ; %%%% calcite burial
% delta_mccb = d13c_A + d_mccb ;
delta_mccb = d13c_A ;

%%%% marine organic
capdelB_0 = 33 ;
Jparam =  5 ;
e_B_co2 = -9 / sqrt(RCO2) ;
e_o2 = Jparam * ( (O/pars.O0) - 1) ;
capdelB = capdelB_0 + e_B_co2 + e_o2 ;
%%%% final calc
d_mocb = d_mccb - capdelB ;
delta_mocb = d13c_A + d_mocb ;

%%%% land plant
capdelP_0 = 19 ;
capdelP = capdelP_0 + e_o2;
%%%% atmospheric
delta_a = (atfrac-1)*( (9483/TEMP_gast)  - 23.89 ) ; 
%%%% final calc
d_locb = delta_a - capdelP ;
delta_locb = d13c_A + d_locb ;

%%%%% S isotopes (copse)
capdelS = 40 ;
delta_mpsb = d34s_S - capdelS ;
delta_evap = 15 ;
delta_pulse_pyr = -30 ;

%%% deltaORG_C*ORG_C 
dy(12) =  locb*(  delta_locb ) + mocb*( delta_mocb )  -   oxidw*delta_G  -   ocdeg*delta_G - Gsub*delta_G  ;

%%% deltaCARB_C*CARB_C 
dy(13) =  mccb*delta_mccb + sfw*delta_mccb  -  carbw*delta_C  - ccdeg*delta_C - Csub*delta_C ;

%%% deltaPYR_S*PYR_S (young)
dy(14) =  mpsb*( delta_mpsb )  - pyrw*delta_PYR  - pyrdeg*delta_PYR - PYRsub*delta_PYR ;

%%% deltaGYP_S*GYP_S (young)
dy(15) =  mgsb*d34s_S   - gypw*delta_GYP  - gypdeg*delta_GYP - GYPsub*delta_GYP ;

%%% delta_A * A
dy(16) = -locb*(  delta_locb ) -mocb*( delta_mocb ) + oxidw*delta_G + ocdeg*delta_G + ccdeg*delta_C + carbw*delta_C - mccb*delta_mccb - sfw*delta_mccb + DOC_ox*-30 + reductant_input*-5 ;

%%% delta_S * S
dy(17) = gypw*delta_GYP + pyrw*delta_PYR -mgsb*d34s_S - mpsb*( delta_mpsb ) + gypdeg*delta_GYP + pyrdeg*delta_PYR + evapdis*delta_evap + pulse_pyr*delta_pulse_pyr ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Strontium system   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% fluxes
Sr_granw = pars.k_Sr_granw *( granw / pars.k_granw ) ;
Sr_basw = pars.k_Sr_basw *( basw / pars.k_basw ) ;
Sr_sedw = pars.k_Sr_sedw *( carbw / pars.k_carbw ) * (SSr/pars.SSr0) ;
Sr_mantle = pars.k_Sr_mantle * DEGASS ;
Sr_sfw = pars.k_Sr_sfw * (sfw/pars.k_sfw) * ( OSr/pars.OSr0 ) ;
Sr_metam = pars.k_Sr_metam * DEGASS * (SSr/pars.SSr0) ;
Sr_sedb = pars.k_Sr_sedb * ( mccb/pars.k_mccb ) * ( OSr/pars.OSr0 ) ;

%%%% fractionation calculations
delta_OSr = y(19) / y(18) ;
delta_SSr = y(21) / y(20) ;

% %%%% original frac
RbSr_bas = 0.1 ;
RbSr_gran = 0.26 ;
RbSr_mantle = 0.066 ;
RbSr_carbonate = 0.5 ;

%%%% frac calcs
dSr0 = 0.69898 ;
tforwards = 4.5e9 + t ;
lambda = 1.4e-11 ;
dSr_bas = dSr0 + RbSr_bas*( 1 - exp(-1*lambda*tforwards) ) ;
dSr_gran = dSr0 + RbSr_gran*( 1 - exp(-1*lambda*tforwards) ) ;
dSr_mantle = dSr0 + RbSr_mantle*( 1 - exp(-1*lambda*tforwards) ) ;

%%%% Ocean [Sr]
dy(18) = Sr_granw + Sr_basw + Sr_sedw + Sr_mantle - Sr_sedb - Sr_sfw ;

%%%% Ocean [Sr]*87/86Sr
dy(19) = Sr_granw*dSr_gran + Sr_basw*dSr_bas + Sr_sedw*delta_SSr + Sr_mantle*dSr_mantle - Sr_sedb*delta_OSr - Sr_sfw*delta_OSr ;

%%%% Sediment [Sr]
dy(20) = Sr_sedb - Sr_sedw - Sr_metam ;

%%%% Sediment [Sr]*87/86Sr
dy(21) = Sr_sedb*delta_OSr - Sr_sedw*delta_SSr - Sr_metam*delta_SSr + SSr*lambda*RbSr_carbonate*exp(lambda*tforwards)  ;

%%%% DOC reservoir
dy(22) = - DOC_ox ;

%%%% uranium system (UNFINISHED)
U_riv = 4.79e7 * (silw/pars.k_silw) ;
U_hydro = 5.7e6 * DEGASS * (U/pars.U0) ;
U_anox = 6.2e6 * (ANOX/0.0025) * (U/pars.U0) ;
U_sed = 3.6e7 * (U/pars.U0) ;

%%%% isotope fractionations
d238U_riv = -0.29 ;
capU_anox = -0.5 ;

%%%% marine U
dy(23) = U_riv - U_hydro - U_anox - U_sed ;
%%%% U * d238U
dy(24) = U_riv*d238U_riv - U_hydro*d238U_sw - U_anox*(d238U_sw - capU_anox) - U_sed*d238U_sw ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Mass conservation check   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_C = A + G + C ;
res_S = S + PYR + GYP ;
iso_res_C = A*d13c_A + G*delta_G + C*delta_C ;
iso_res_S = S*d34s_S + PYR*delta_PYR + GYP*delta_GYP ;


%%%%% istopic composition of inputs
d_in = ( oxidw*delta_G + ocdeg*delta_G + ccdeg*delta_C + carbw*delta_C + DOC_ox*-30 + reductant_input*-5 ) / ( oxidw + ocdeg + ccdeg + carbw + DOC_ox + reductant_input )  ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Print full states for single run   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sensanal == 0
    workingstate.iso_res_C(stepnumber,1) = iso_res_C ;
    workingstate.iso_res_S(stepnumber,1) = iso_res_S ;
    workingstate.res_C(stepnumber,1) = res_C ;
    workingstate.res_S(stepnumber,1) = res_S ;
    workingstate.time(stepnumber,1) = t;
    workingstate.temperature(stepnumber,1) = TEMP_gast ;
    workingstate.tempC(stepnumber,1) = TEMP_gast - 273 ;
    workingstate.P(stepnumber,1) = P ;
    workingstate.O(stepnumber,1) = O ;
    workingstate.A(stepnumber,1) = A ;
    workingstate.S(stepnumber,1) = S ;
    workingstate.G(stepnumber,1) = G ;
    workingstate.C(stepnumber,1) = C ;
    workingstate.PYR(stepnumber,1) = PYR ;
    workingstate.GYP(stepnumber,1) = GYP ;
    workingstate.CAL(stepnumber,1) = CAL ;
    workingstate.N(stepnumber,1) = N ;
    workingstate.OSr(stepnumber,1) = OSr ;
    workingstate.SSr(stepnumber,1) = SSr ;
    %%%%%%% print isotope information
    workingstate.d13c_A(stepnumber,1) = d13c_A ;
    workingstate.delta_mccb(stepnumber,1) = delta_mccb ;
    workingstate.d34s_S(stepnumber,1) = d34s_S ;
    workingstate.delta_G(stepnumber,1) = delta_G ;
    workingstate.delta_C(stepnumber,1) = delta_C ;
    workingstate.delta_PYR(stepnumber,1) = delta_PYR ;
    workingstate.delta_GYP(stepnumber,1) = delta_GYP ;
    workingstate.delta_o(stepnumber,1) = delta_o ;
    workingstate.delta_a(stepnumber,1) = delta_a ;
    workingstate.d_mccb(stepnumber,1) = d_mccb ;
    workingstate.d_locb(stepnumber,1) = d_locb;
    workingstate.d_mocb(stepnumber,1) = d_mocb ;
    workingstate.delta_OSr(stepnumber,1) = delta_OSr ;
    %%%%%%% print forcings
    workingstate.UPLIFT(stepnumber,1) = UPLIFT ;
    workingstate.DEGASS(stepnumber,1) = DEGASS ;
    workingstate.W(stepnumber,1) = W ;
    workingstate.EVO(stepnumber,1) = EVO ;
    workingstate.CPland(stepnumber,1) = CPland_relative ;
    workingstate.Bforcing(stepnumber,1) = Bforcing ;
    % workingstate.SOLAR(stepnumber,1) = SOLAR ;
    workingstate.BAS_AREA(stepnumber,1) = BAS_AREA ;
    workingstate.GRAN_AREA(stepnumber,1) = GRAN_AREA ;
    workingstate.CARB_AREA(stepnumber,1) = CARB_AREA ;
    workingstate.PG(stepnumber,1) = PG ;
    %%%%%%%% print variables
    workingstate.RCO2(stepnumber,1) = RCO2 ;
    workingstate.RO2(stepnumber,1) = RO2 ;
    workingstate.mrO2(stepnumber,1) = mrO2 ;
    workingstate.VEG(stepnumber,1) = VEG ;
    workingstate.ANOX(stepnumber,1) = ANOX ;
    % workingstate.ALBEDO(stepnumber,1) = ALBEDO ;
    %%%%%%%% print fluxes
    workingstate.mocb(stepnumber,1) = mocb ;
    workingstate.locb(stepnumber,1) = locb ;
    workingstate.mccb(stepnumber,1) = mccb ;
    workingstate.mpsb(stepnumber,1) = mpsb ;
    workingstate.mgsb(stepnumber,1) = mgsb ;
    workingstate.silw(stepnumber,1) = silw ;
    workingstate.carbw(stepnumber,1) = carbw ;
    workingstate.oxidw(stepnumber,1) = oxidw ;
    workingstate.basw(stepnumber,1) = basw ;
    workingstate.granw(stepnumber,1) = granw ;
    workingstate.phosw(stepnumber,1) = phosw ;
    workingstate.psea(stepnumber,1) = psea ;
    workingstate.nfix(stepnumber,1) = nfix ;
    workingstate.denit(stepnumber,1) = denit ;
    workingstate.VEG(stepnumber,1) = VEG ;
    % workingstate.phosw_bas(stepnumber,1) = phosw_bas; 
    % workingstate.phosw_gran(stepnumber,1) = phosw_gran; 
    % workingstate.phosw_carb(stepnumber,1) = phosw_carb; 
    % workingstate.phosw_oxid(stepnumber,1) = phosw_oxid; 
    workingstate.pyrw(stepnumber,1) = pyrw ;
    workingstate.gypw(stepnumber,1) = gypw ;
    workingstate.ocdeg(stepnumber,1) = ocdeg ;
    workingstate.ccdeg(stepnumber,1) = ccdeg ;
    workingstate.pyrdeg(stepnumber,1) = pyrdeg ;
    workingstate.gypdeg(stepnumber,1) = gypdeg ;
    workingstate.sfw(stepnumber,1) = sfw ;
    workingstate.DOC_res(stepnumber,1) = DOC_res ;
    workingstate.evapdis(stepnumber,1) = evapdis ;
    workingstate.pulse_pyr(stepnumber,1) = pulse_pyr ;
    workingstate.DOC_ox(stepnumber,1) = DOC_ox ;
    workingstate.Sr_granw(stepnumber,1) = Sr_granw ;
    workingstate.Sr_basw(stepnumber,1) = Sr_basw ;
    workingstate.Sr_sedw(stepnumber,1) = Sr_sedw ;
    workingstate.Sr_mantle(stepnumber,1) = Sr_mantle ;
    workingstate.dSSr(stepnumber,1) = dSSr ;
    workingstate.relativenewp(stepnumber,1) = newp/pars.newp0 ;
    workingstate.U(stepnumber,1) = U ;
    workingstate.d238U_sw(stepnumber,1) = d238U_sw ;
    workingstate.d_in(stepnumber,1) = d_in ;
    %%%%%%%% print time
    workingstate.time_myr(stepnumber,1) = t_geol ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Print plotting states only in sensanal   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sensanal == 1
    workingstate.delta_mccb(stepnumber,1) = delta_mccb ;
    workingstate.d34s_S(stepnumber,1) = d34s_S ;
    workingstate.delta_OSr(stepnumber,1) = delta_OSr ;
    workingstate.SmM(stepnumber,1) = 28*S/pars.S0 ;
    workingstate.CO2ppm(stepnumber,1) = RCO2*280 ;
    workingstate.mrO2(stepnumber,1) = mrO2 ;
    workingstate.T_gast(stepnumber,1) = TEMP_gast - 273 ;
    workingstate.time_myr(stepnumber,1) = t_geol ;
    workingstate.time(stepnumber,1) = t;
    workingstate.ANOX(stepnumber,1) = ANOX ;
    workingstate.DOC_res(stepnumber,1) = DOC_res ;
end

%%%%%%% output timestep if specified
if pars.telltime ==1
    if mod(stepnumber,pars.display_resolution) == 0 
        %%%% print model state to screen
        fprintf('Model step: %d \t', stepnumber); fprintf('time: %d \n', t_geol)
    end
end



%%%% final action record current model step
stepnumber = stepnumber + 1 ;



%%%% option to bail out if model is running aground
if stepnumber > pars.bailnumber
   terminateExecution
end




end



