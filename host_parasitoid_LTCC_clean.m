% Matlab code to simulate (stochastic) delay differential equation, stage-structured, 
% Host-Parasitoid models under long-term environmental change.
%
% Temperature drives changes in Host & Parasitoid larval development time,
% reduce Parasitoid attack rate and increase Parasitoid larval mortality 
% rate, all based on empirical measurements.
%
% Code originally written by Christophe Coste <christophe.coste@nina.no>, 
% <c.f.d.coste.ecology@gmail.com>, funded by NERC grant number
% NE/W006731/1 awarded to Mike S Fowler, Steve M Sait and Chenggui Yuan.
% Developed from a previous stochastic HP model version by David Gilljam
% <david.gilljam@slu.se>, funded by NERC grant number NE/N002849/1 awarded 
% to Mike S Fowler, Steve M Sait and Chenggui Yuan.
% This version cleaned up by by Mike S Fowler <m.s.fowler@swansea.ac.uk> 
% 03 Nov 2026.
% Runs successfully on Matlab R2023b (v23.2) and R2025a on 64-bit Mac Apple 
% M3 Pro. Doesn't require additional Toolboxes.

% N.B. Fix required: Parasitoid Adult densities should be reset to 0 after
% declining below ymin (e.g., ymin = 1e-6), but remain at positive value
% just below this pseudo-extinction threshold.

clear; clc; close all;

%% ===================== Main settings =====================
global d_HL a_ref a_min Tref_attack I_HAsize I_HA_START I_PAsize I_PA_START P_threshold

moving_average_only = false;    % true: larval durations from the moving average only (one model run, faster)
                                % false: also compute them with the dde, plot both (L158, currently commented out), and run the model with each

time = 550;  % total simulation time, 550 matches 1.8 year experimental time for Li et al (in prep)

% temperature time series: long-term temperature increase (alternative scenarios commented)
T = linspace(28,30,time);

d_HL = 0;   % Background host larval mortality (per day), in absence of parasitoids

% parasitoid attack rate: a(T) = max(0, a_ref + a_Tslope*(T-Tref_attack)), falling from a_ref to a_min over the simulation run
a_ref       = 0.015;          % Attack rate at reference temperature (mean of small & large larval attack rates from Wearing et al (2004, Am Nat), or White et al (2007) = 0.01
a_min       = 0.5.*a_ref;     % Min Par. attack rate at the end of the experiment.
Tref_attack = 28;             % Reference temperature (deg C)

% parasitoid mortality: d(T) = max(0, d_ref + d_Tslope*(T - Tref))
d_PL_ref    = 0.1;    % Parasitoid juvenile mortality at reference temperature (per day)
d_PL_Tslope = 0.05;   % Change in juvenile mortality per deg C
Tref_dPL    = 28;     % Reference temperature (deg C)
d_PA_ref    = 0.1;    % Parasitoid adult mortality at reference temperature (per day)
d_PA_Tslope = 0;      % Change in adult mortality per deg C
Tref_dPA    = 28;     % Reference temperature (deg C)

% introductions
I_HAsize   = 30;   % Starting number of Hosts
I_HA_START = 0;    % Starting time for introduction of Hosts
I_PAsize   = 3;    % Starting number of Parasitoids
I_PA_START = 25;   % Starting time for introduction of parasitoids (rough average of day 25 vs 146)

% parasitoid extinction
P_threshold = 1e-6;   % parasitoids (adults + larvae) below this are considered extinct

%% parameters Host larval growth rate as function of temperature
% Coefficients estimated from single-generation experiments by Mugabo et al (unpublished)
global T_a T_b1 T_b2
T_a = -0.279045;
T_b1 = 0.023349;
T_b2 = -4.02e-04;

alpha = -T_b2;
Topt = T_b1/(2*alpha);      % Optimal development temperature for Host

%% parameters Parasitoid larval growth rate as function of temperature
% Coefficients estimated from single-generation experiments by Mugabo et al (unpublished)
global T_aPL T_b1PL T_b2PL
T_aPL = -0.2336379489;
T_b1PL = 0.019642587;
T_b2PL = -3.425152e-04;

alphaPL = -T_b2PL;
ToptPL = T_b1PL/(2*alphaPL);  % Optimal development temperature for Parasitoid

%%  growth rate time series
g=gtemp(T);
gPL=gtempPL(T);

Tsto=griddedInterpolant(1/2:1:(time-1/2),T,'nearest');

global tempsto gHLoft gPLoft
tempsto=@(t) Tsto(t);
gHLoft=@(t) gtemp(tempsto(t));
gPLoft=@(t) gtempPL(tempsto(t));

%%  larvae stage duration time series

global tauinit taufinal tauPLinit tauPLfinal
tauinit=   1/g(1);
taufinal=1/g(end);
tauPLinit=   1/gPL(1);
tauPLfinal=1/gPL(end);

%% Generate the larval developmental delay times, tau proxy via "moving average" 
% approximation (Coste et al, in prep), with a time step of 1/dilat day

dilat=100;
gdil=repmat(g,dilat,1);
gdil=gdil(:)';
gPLdil=repmat(gPL,dilat,1);
gPLdil=gPLdil(:)';

tt1=round(tauinit*dilat); tt=tt1;  %initial window
ttPL1=round(tauPLinit*dilat);ttPL=ttPL1; %initial window

for t=tt1:time*dilat
    mgma(t)=(sum(gdil(t-(tt-1):t)))/tt;
    tt=round(dilat/mgma(t));
end
for t=1:tt1-1
    mgma(t)=((tt1-t)*gdil(1) + sum(gdil(1:t)))/tt1;
end

for t=ttPL1:time*dilat
    mgmaPL(t)=(sum(gPLdil(t-(ttPL-1):t)))/ttPL;
    ttPL=round(dilat/mgmaPL(t));
end
for t=1:ttPL1-1
    mgmaPL(t)=((ttPL1-t)*gPLdil(1) + sum(gPLdil(1:t)))/ttPL1;
end

invmgmadilat=1./mgma;
invmgmaPLdilat=1./mgmaPL;

invmgma2b=griddedInterpolant(1:1:time*dilat,invmgmadilat);
invmgma2PLb=griddedInterpolant(1:1:time*dilat,invmgmaPLdilat);

global invmgmaoft2 invmgmaPLoft2
invmgmaoft2=@(t) invmgma2b(t*dilat);
invmgmaPLoft2=@(t) invmgma2PLb(t*dilat);


%% Set dde tolerances

optsstrong = ddeset('RelTol',0.4e-07,'AbsTol',1e-10); %minimum
optsdavid = ddeset('RelTol',1e-05,'AbsTol',1e-08); %Davids parameters

opts = optsstrong;   % used for the dde larval durations below

timevec=0:time;

%% tau computed via sdde (only if moving_average_only = false)
if ~moving_average_only
    tic
        soltau=ddesd(@tauHLt2,@taudelays,tauinit,[-100, time],opts);
    toc

    tic
        soltauPL=ddesd(@tauPLt2,@taudelays,tauPLinit,[-100, time],opts);
    toc

    taut=deval(soltau,timevec);
    tauPLt=deval(soltauPL,timevec);
 
  %   % compare the larval durations from the two methods
  %   figure(10)
  %   subplot(211)
  %       plot(timevec,taut,'b'); hold on
  %       plot(timevec,invmgmaoft2(timevec),'g--')
  %       ylabel('Duration (days)')
  %       title('host larval duration')
  %       legend('dde','moving average')
  %   subplot(212)
  %       plot(timevec,tauPLt,'b'); hold on
  %       plot(timevec,invmgmaPLoft2(timevec),'g--')
  %       xlabel('Time (days)'); ylabel('Duration (days)')
  %       title('parasitoid larval duration')
  %       legend('dde','moving average')
  % fprintf('Larval development time: largest difference between the dde and the moving average over the run: host %.2g days, parasitoid %.2g days\n', ...
  %       max(abs(taut - invmgmaoft2(timevec))), max(abs(tauPLt - invmgmaPLoft2(timevec))))
end

%% larval durations used by the model (set in the run loop below)
global tauoft tauPLoft

%% Host Alone System
% parameters: fixed duration, birth rate, death rates, and density dependence
global tau_HE    tau_HP    tau_HA  b c d_HE d_HP d_HA
tau_HE    = 4.3;      % Duration of host egg stage
tau_HP    = 7;        % Duration of host puppae stage
tau_HA    = 5.5;      % Duration of host adult stage
b         = 35;       % Daily host adult fecundity (no. eggs/day)
c         = 0.00007;          % Competition based on H-alone parameterisation.
d_HE        = 0.017;   % Background host egg mortality (per day)
d_HP        = 0.1;       % Background host PUPAL mortality
d_HA        = 0.1;       % Background host adult mortality

% computation of stage survival rates from the input rates
global sigma_HE sigma_HPup sigma_HA
sigma_HE    = exp(-d_HE * tau_HE) ;   % Time independent survival probabilities
sigma_HPup  = exp(-d_HP * tau_HP);
sigma_HA    = exp(-d_HA * tau_HA);

%% with Parasitoid

global tau_PA  k a_Tslope
tau_PA    = 1;    % Duration of parasitoid adult stage; MSF: updated to 1 (was 2) following advice from SMS
k         = 1;    % Interference parameter; MSF: updated to 1 (was 0.01) to match Wearing et al (2004, Am Nat)

% slope of the attack rate (per deg C)
if (var(T)==0)
    a_Tslope = 0;     % If there is no temperature change, revert to constant attack rate.
else
    a_Tslope = -(a_ref-a_min)/(T(end)-T(1));   % MSF: Wearing et al's (2004) -[large-small larval attack rates]/XºC change
end

% Temperature-dependent parasitoid juvenile (larval) mortality
global dPLoft DPLcum
dPLoft = @(t) max(0, d_PL_ref + d_PL_Tslope*(tempsto(t) - Tref_dPL));

% Cumulative juvenile mortality, for survival through the parasitoid larval stage
% (temperature is constant within each day, so this integral is exact)
dPLdaily    = max(0, d_PL_ref + d_PL_Tslope*(T - Tref_dPL));
DPLcum_grid = griddedInterpolant(0:time, [0 cumsum(dPLdaily)], 'linear', 'linear');
DPLcum      = @(t) DPLcum_grid(t);

% Temperature-dependent parasitoid adult mortality
global dPAoft DPAcum
dPAoft = @(t) max(0, d_PA_ref + d_PA_Tslope*(tempsto(t) - Tref_dPA));

% Cumulative adult mortality, for survival through the adult stage
% (temperature is constant within each day, so this integral is exact)
dPAdaily    = max(0, d_PA_ref + d_PA_Tslope*(T - Tref_dPA));
DPAcum_grid = griddedInterpolant(0:time, [0 cumsum(dPAdaily)], 'linear', 'linear');
DPAcum      = @(t) DPAcum_grid(t);


%% model with Parasitoid Larvae implemented
% Run once with the moving-average larval durations, and a second time with
% the dde durations if moving_average_only = false
if moving_average_only
    tau_methods = {'moving average'};
else
    tau_methods = {'moving average', 'dde'};
end

global P_extinct
opts_ext = ddeset(optsdavid, 'Events', @parasitoid_extinct);

% Break the run where the parasitoid introduction switches on/off, so the solver
% does not step across these jumps (stepping across them gave small negative values)
t_breaks = [0, I_PA_START, I_PA_START+1, I_PA_START+1+tau_PA];

for m = 1:numel(tau_methods)
    if strcmp(tau_methods{m}, 'moving average')
        tauoft = invmgmaoft2;   tauPLoft = invmgmaPLoft2;                      % moving-average durations
    else
        tauoft = @(t) deval(soltau,t);   tauPLoft = @(t) deval(soltauPL,t);    % dde durations
    end
    P_extinct = false;

    tic
    sol = ddesd(@ddesdfullsystem_wPL,@ddesdfullsystemdelays,@ddesdfullsystemhistory_wPL,[t_breaks(1), t_breaks(2)],optsdavid);
    for j = 2:numel(t_breaks)-1
        sol = ddesd(@ddesdfullsystem_wPL,@ddesdfullsystemdelays,sol,[t_breaks(j), t_breaks(j+1)],optsdavid);
    end
    sol = ddesd(@ddesdfullsystem_wPL,@ddesdfullsystemdelays,sol,[t_breaks(end), time],opts_ext);
    t_extinct = NaN;
    if sol.x(end) < time                % parasitoids fell below the threshold before the end
        t_extinct = sol.x(end);          % day of parasitoid extinction
        P_extinct = true;               % from now on parasitoid abundances stay at 0
        y0 = sol.y(:,end);
        y0([4 5]) = 0;                  % parasitoid adults and larvae set to 0
        sol = ddesd(@ddesdfullsystem_wPL,@ddesdfullsystemdelays,sol,[sol.x(end), time],ddeset(optsdavid,'InitialY',y0));
    end
    solsystem_wPL = sol;
    toc

    solsystemt_wPL = deval(solsystem_wPL,0:time);
    % rows: 1 host larvae, 2 host adults, 3 host larval survival, 4 parasitoid adults, 5 parasitoid larvae

    %% adult host and adult para
    figure(m);
    subplot(311);   % Adult H & P abundances, log scale
    semilogy(solsystemt_wPL([2,4],:)');
    ylabel('Adult Abundances')
    title('abundance time series')
    legend('host adults','para adults')

    subplot(312);   % Adult H & P abundances, linear scale
    plot(solsystemt_wPL([2,4],:)');
    ylabel('Adult Abundances')
    title('abundance time series')
    legend('host adults','para adults')

    subplot(313);   % Larval H & P abundances, log scale
    semilogy(solsystemt_wPL([1,5],:)');
    xlabel('Time (days)');
    ylabel('Larval Abundances')
    legend('host larvae','para larvae')

    sgtitle(['Larval durations: ' tau_methods{m}])
    results{m}   = solsystemt_wPL;   % output of each run (1 = moving average, 2 = dde)
    t_ext_all(m) = t_extinct;        % extinction day of each run (NaN = no extinction)
end


%% Subfunctions g: empirically derived Host & Parasitoid larval growth rates

function gtempout =gtemp(temp)
global T_a T_b1 T_b2
   gtempout =  T_a+  T_b1*temp +  T_b2*(temp.^2);
   gtempout(gtempout < 0) = 0; % Don't allow negative growth-rates
end

function gtempout =gtempPL(temp)
global T_aPL T_b1PL T_b2PL
   gtempout =  T_aPL+  T_b1PL*temp +  T_b2PL*(temp.^2);
   gtempout(gtempout < 0) = 0; % Don't allow negative growth-rates
end


%% functions tau standalone (dde for the larval durations)

function dydt = tauHLt2(t,y,Z)
    global gHLoft
    dydt=1-(gHLoft(t)/gHLoft(t-y(1)));
end

function dydt = tauPLt2(t,y,Z)
    global gPLoft
    dydt=1-(gPLoft(t)/gPLoft(t-y(1)));
end

function s = taudelays(t,y)
    s=t-0;
end


%% functions host-parasitoid system

function s = ddesdfullsystemhistory_wPL(t)
    global c k a_ref tauinit d_HL
    HL=0;
    HA=0;
    PA=0;
    HLsurvival= exp(-tauinit*( (c*HL + d_HL) + k*log(1+(a_ref/k)*PA)));
    PL=0;
    s=[HL, HA, HLsurvival, PA, PL]';
end


function dydt = ddesdfullsystem_wPL(t,y,Z)
    global   tau_HP    tau_HA  b c d_HA  d_HL  sigma_HE  sigma_HPup    sigma_HA     ;
    global   k a_ref a_Tslope Tref_attack tau_PA ;
    global gHLoft tauoft;
    global  gPLoft tauPLoft tempsto;
    global dPAoft DPAcum;
    global dPLoft DPLcum;
    global P_extinct;
    
    tauHLoft = tauoft;
    
    HL   = 1;    % Host larvae
    HA   = 2;    % Host adults
    SHL  = 3;    % Host larvae dens.dep.prob survival
    PA   = 4;    % Parasitoid adults
    PL   = 5;    % Parasitoid larvae
    
    % Temperature-dependent attack rates evaluated at the biological time of attack
    a_now = max(0, a_ref + a_Tslope*(tempsto(t)-Tref_attack));
    t_HLdelay = t - tauoft(t);
    t_PLdelay1 = t - tauPLoft(t);
    t_PLdelay2 = t - tau_PA - tauPLoft(t-tau_PA);
    
    % For times before the simulated temperature series starts, use the reference rate.
    % (This is consistent with the constant-temperature prehistory.)
    a_HLdelay = a_ref;
    if t_HLdelay >= 0.5
        a_HLdelay = max(0, a_ref + a_Tslope*(tempsto(t_HLdelay)-Tref_attack));
    end
    a_PLdelay1 = a_ref;
    if t_PLdelay1 >= 0.5
        a_PLdelay1 = max(0, a_ref + a_Tslope*(tempsto(t_PLdelay1)-Tref_attack));
    end
    a_PLdelay2 = a_ref;
    if t_PLdelay2 >= 0.5
        a_PLdelay2 = max(0, a_ref + a_Tslope*(tempsto(t_PLdelay2)-Tref_attack));
    end
    
    % death rate of host larvae: competition + background mortality + parasitoids
    D_HL_T =  (c*y(HL) + d_HL) + k*log(1+(a_now/k)*y(PA)) ;   % c*Nl(t) + d_HL + k.ln(1+(a.Nap(t))/k)
    D_HL_Tdelay = (c*Z(HL, 7) + d_HL) + k*log(1+(a_HLdelay/k)*Z(PA,7)); % same but computed with abundances with delay tauHL
    Para_Tdelay1 = k*log(1+(a_PLdelay1/k)*Z(PA,8));
    Para_Tdelay2 = k*log(1+(a_PLdelay2/k)*Z(PA,9));
    D_HL_T1 = k*log(1+(a_now/k)*y(PA)) ;
    
    % Host larvae
    R_HL  = b * Z(HA, 1) * sigma_HE ; %recruitment :bNa(t-τ_e) e^(-de.τe)
    M_HL  = b * Z(HA, 2) * sigma_HE * y(SHL) * gHLoft(t)/gHLoft(t-tauHLoft(t)); %maturation : bNa(t-tau(t)-τ_e)e^(-de.τe)Sl(t)(gL(t))/(gL(t-tau(t)) )
    D_HL  = y(HL) * D_HL_T;  %instantaneous death rate of host larvae:
    dHLdt = R_HL - M_HL - D_HL;
    
    % Host adults
    R_HA  = b * Z(HA, 3) * sigma_HE * Z(SHL, 4) * gHLoft(t-tau_HP)/gHLoft(t-tau_HP-tauHLoft(t-tau_HP)) * sigma_HPup + I_HA(t);
    % recruitment of host adults:bNa(t-τpp-tau(t-τpp)-τe)e^(-deτe)Sl(t-τpp)(gL(t-τpp))/(gL(t-τpp-tau(t-τpp)))e^(-dppτpp)
    M_HA  = (b * Z(HA, 5) * sigma_HE * Z(SHL, 6) * gHLoft(t-tau_HA-tau_HP)/gHLoft(t-tau_HA-tau_HP-tauHLoft(t-tau_HA-tau_HP))  * sigma_HPup+ I_HA(t-tau_HA)) * sigma_HA;
    %maturation of host adults bNa(t-τa-τpp-tau(t-τa-τpp)-τe) e^(-deτe)Sl(t-τa-τpp)(gL(t-τa-τpp))/(gL(t-τa-τpp-tau(t-τa-τpp)))e^(-dppτpp)e^(-daτa)
    D_HA  = d_HA * y(HA);
    dHAdt = R_HA - M_HA - D_HA;
    
    % Through-stage survival of larvae (Sl)
    dSHLdt = y(SHL) * ( D_HL_Tdelay * gHLoft(t)/gHLoft(t-tauHLoft(t)) - D_HL_T );
    
    % Parasitoid adults
    sigma_PL_now = exp(-(DPLcum(t) - DPLcum(t-tauPLoft(t))));   % survival through the larval stage for larvae maturing now (entered at t-tauPLoft(t))
    R_PA  = sigma_PL_now*Para_Tdelay1*Z(HL,8) * gPLoft(t)/gPLoft(t-tauPLoft(t)) + I_PA(t);
    %recruitment : (survival through larval stage) * k.ln(1+(a.Nap(t-τauP(t)))/k) * Nl(t-τauP(t)) * (gP(t))/(gP(t-τauP(t)) )
    
    sigma_PA_t = exp(-(DPAcum(t) - DPAcum(t-tau_PA)));   % survival through the adult stage for adults leaving now (entered at t-tau_PA)
    t_PAin = t - tau_PA;   % time at which the adults leaving now became adults
    sigma_PL_old = exp(-(DPLcum(t_PAin) - DPLcum(t_PAin-tauPLoft(t_PAin))));   % their survival through the larval stage
    M_PA  =  (sigma_PL_old * Para_Tdelay2 * Z(HL,9) * gPLoft(t-tau_PA)/gPLoft(t-tau_PA-tauPLoft(t-tau_PA))+  I_PA(t-tau_PA) )*sigma_PA_t   ;
    %maturation is same than recruitment observed at t-tau_PA  * thru-stage survival for PA
    D_PA  = dPAoft(t) * y(PA);
    dPAdt = R_PA - M_PA - D_PA;
    
    % Parasitoid larvae
    R_PL = D_HL_T1*y(HL);
    M_PL = R_PA;
    D_PL = dPLoft(t) * y(PL);
    dPLdt = R_PL - M_PL - D_PL;
    
    dydt = [dHLdt dHAdt dSHLdt dPAdt dPLdt ]';
    
    if P_extinct
        dydt([PA PL]) = 0;   % parasitoids extinct: adults and larvae stay at 0
    end
    
    % abundances below ymin are not allowed to decline further
    ymin = 10^-6;
    dydt((y < ymin) & (dydt < 0)) = 0;

end

function d = ddesdfullsystemdelays(t,y)
global tau_HE tau_HP tau_HA ;
global  tauoft;
global tau_PA;
global  tauPLoft;

  l1  = t - tau_HE;
  l2  = t - tau_HE - tauoft(t);
  l3  = t - tau_HP-tau_HE-tauoft(t-tau_HP);
  l4  = t - tau_HP;
  l5  = t - tau_HP-tau_HE-tau_HA-tauoft(t-tau_HP-tau_HA);     %  t - tau_HE - tau_HL(t - tau_HP-tau_HA) - tau_HP - tau_HA
  l6  = t - tau_HP - tau_HA;
  l7 = t  - tauoft(t);
  % parasitoid
  l8= t - tauPLoft(t);
  l9 = t - tau_PA - tauPLoft(t-tau_PA);

  d = [l1 l2 l3 l4 l5 l6 l7 l8 l9]';
end


%% Subfunction inoculations/initialising Host/Parasitoid populations

function y=I_HA(t)
global I_HAsize I_HA_START
I_HA_STOP =I_HA_START+ 1;
y=(I_HA_START<=t).*(t<I_HA_STOP).*I_HAsize ;
end

function y=I_PA(t)
global I_PAsize I_PA_START
I_PA_STOP =I_PA_START +1 ;
y=(I_PA_START<=t).*(t<I_PA_STOP).*I_PAsize ;
end


%% Event handling: parasitoid extinction

function [value, isterminal, direction] = parasitoid_extinct(t, y, Z)
global P_threshold
value      = y(4) + y(5) - P_threshold;   % total parasitoids minus threshold
isterminal = 1;                           % stop the integration
direction  = -1;                          % only when falling below the threshold
end
