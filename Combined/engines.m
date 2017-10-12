clc; close all; clear all;

%% %%%%% 2 %%%%%%%
y = 1.4
R = 287; %j/kg-K
cp = y*R/(y-1);
% a -> diff. -> 1 -> fan -> 2 -> HP_comp -> 3 -> Comb -> 4 -> HP_turb -> 5 
%     -> LP_turb -> 6 -> Nozzel -> 7

M = 0.85;
pa = 23800% Pa
rhoa = 0.38% kg/m3
Ta = 218.8% K

HP_comp_r = 20
while HP_comp_r < 61
    
BPR = linspace(5,15,11)


%Atmosphere
p0a = pa*(1+(y-1)*M*M/2)^(y/(y-1));
T0a = Ta*(1+(y-1)*M*M/2);

%Diffuser
p01 = p0a;
T01 = T0a;

%fan
fan_r = 1.5;
p02 = fan_r*p01;
n_fan = 0.9;
T02 = T01*(1+(1/n_fan)*(fan_r^((y-1)/y)-1));

%HP compressor
%HP_comp_r=[20 60]
p03 = HP_comp_r*p02;
n_comp = 0.9;
T03 = T02*(1+(1/n_comp)*(HP_comp_r^((y-1)/y)-1));

%Combustor
% phi = 0.25;
LHV = 43e6; %j/kg
% C12H26 + 18.5 (O2 + 3.76N2) -> 12 CO2 + 13 H2O + 69.56 N2
FoA = 170.34/(18.5*32+18.5*3.76*28);% fuel to air ratio
p04 = p03;
% T04 = T03 + phi*FoA*LHV/cp;

T04 = 1600;
phi = (T04-T03)*cp/(FoA*LHV)
%HP turbine
n_turb = 0.95;
T05 = T04 - (T02/n_comp)*(((HP_comp_r)^((y-1)/y))-1);
p05 = p04*(1-(1/(n_comp*n_turb*T04))*T02*(((HP_comp_r)^((y-1)/y))-1))^(y/(y-1));

%Lp turbine
T06 = T05 - (1+BPR)*(T01/n_fan)*(((fan_r).^((y-1)/y))-1);
p06 = p05*(1-((1+BPR)/(n_fan*n_turb*T05))*T01*(((fan_r).^((y-1)/y))-1)).^(y/(y-1));

%nozzel
Ve_c = sqrt(2*(y/(y-1))*R*T06.*(1-(pa./p06).^((y-1)/y)));
Ve_f = sqrt(2*(y/(y-1))*R*T02.*(1-(pa./p02).^((y-1)/y)));

sos = sqrt(1.4*pa/rhoa);
V = sos*M;
Thrust = (BPR*(Ve_f-V)+(Ve_c-V))
Np = Thrust.*V./(BPR*((Ve_f.^2)/2-(V.^2)/2)+((Ve_c.^2)/2-(V.^2)/2))
TSFC = 3600*4.44*2.2*phi*FoA./Thrust


figure(1)
subplot(1,3,1)
plot(BPR,Thrust); hold on
xlabel('BPR')
ylabel('Specific Thrust (N-s/kg)')
subplot(1,3,2)
plot(BPR,Np); hold on
ylabel('Propulsive Efficiency')
xlabel('BPR')
subplot(1,3,3)
plot(BPR,TSFC); hold on
ylabel('TSFC (lb/hr-lb)')
xlabel('BPR')

% 
HP_comp_r = HP_comp_r + 10

end
lgd = legend('20','30','40','50','60');

% This time around, the compressor ratio reduces the specific thrust
% because it as the compression ratio increases, the specific energy that
% can be added is limited by T04. On the other hand, TSFC still benefits
% from a higher compression ratio as the cycle is thermodynamically more
% efficient at higher pressures. Once again, the compressor ratio has no
% effect on the peak propulsive effiency.
% The HP compressor ratio and BPR selected where 40 and 13 respectively.
% The compressor ratio was selected because it proves to be a good middle
% ground between specifc thrust and TSFC. The BPR was chosen because it is
% the most efficent and 'powerful' setup at that compression ratio.

G_thrust = 105000*4.44822% N
Cruise_thrust = 0.15*G_thrust  % 70059 N
specific_thrust = 1455 % N-s/kg through the core (@ BPR = 13, HP_r = 40
BPR = 13
mflow = (Cruise_thrust/specific_thrust)*(BPR+1)  %674 kg/s

vflow = mflow/rhoa

area = vflow/V

diameter = sqrt(area*4/pi())% 2.997m

% This diameter is still quite similar to that of modern larger BPR
% Turbofans



%%
clc; clear all;
TO_thrust = 105000*4.44822*0.5 % required thrust at takeof
BPR = 13
M = 200/761;
pa = 101350% Pa
rho = 1.225% kg/m3
Ta = 288% K
V_TO = 200*0.44704
TO_area = 7.04
TO_vflow = TO_area*V_TO
TO_mflow = TO_vflow*rho
specific_thrust_required =(TO_thrust/TO_mflow)*(BPR+1) %4240N

%T04 required is about 2050K with phi = 0.36

% A bit hot but takeoff only occurs for a limited time so it doesnt matter
% if the engine melts a little bit for short periods of time (at least
% thats what was hinted during lecture).

