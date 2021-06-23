%%%% make forcing files for COPSE reloaded

%%%% load workspace
load reloaded_forcings.mat

%%%% interpolate all onto standard t grid
tgrid = [-1000 : 1 : 0] ;
new_B = interp1( -1.*t_B , B , tgrid ) ;
new_BA = interp1( -1.*t_BA , BA , tgrid ) ;
new_Ca = interp1( -1.*t_Ca , Ca , tgrid ) ;
new_D = interp1( -1.*t_D , D , tgrid ) ;
new_E = interp1( [ -1000 -1.*t_E] , [ 0 E ] , tgrid ) ;
new_GA = interp1( -1.*t_GA , GA , tgrid ) ;
new_PG = interp1( -1.*t_PG , PG , tgrid ) ;
new_U = interp1( -1.*t_U , U , tgrid ) ;
new_W = interp1( -1.*t_W , W , tgrid ) ;
new_coal = interp1( -1.*t_coal , coal , tgrid ) ;
new_epsilon = interp1( -1.*t_epsilon , epsilon , tgrid ) ;
CP(1:400) = 1 ;
new_CP = interp1(  t_CP  , CP  , tgrid ) ;


%%%% make forcings structure
forcings.t = tgrid ;
forcings.B = new_B ;
forcings.BA = new_BA ;
forcings.Ca = new_Ca ;
forcings.CP = new_CP ;
forcings.D = new_D ;
forcings.E = new_E ;
forcings.GA = new_GA ;
forcings.PG = new_PG ;
forcings.U = new_U ;
forcings.W = new_W ;
forcings.coal = new_coal ;
forcings.epsilon = new_epsilon ;

%%%% plot all
figure
hold on
plot( forcings.t , forcings.B ) 
plot( forcings.t , forcings.BA ) 
plot( forcings.t , forcings.Ca ) 
plot( forcings.t , forcings.CP ) 
plot( forcings.t , forcings.D ) 
plot( forcings.t , forcings.E ) 
plot( forcings.t , forcings.GA ) 
plot( forcings.t , forcings.PG ) 
plot( forcings.t , forcings.U ) 
plot( forcings.t , forcings.W ) 
plot( forcings.t , forcings.coal ) 
plot( forcings.t , forcings.epsilon ) 


