clc; clear all; close all;
M = 3;
pa = 17873.9% Pa
rhoa = 0.287407% kg/m3
Ta = 216.650% K
y = 1.4
R = 287; %j/kg-K
cp = y*R/(y-1);
sos = sqrt(1.4*pa/rhoa);
LHV = 43e6; %j/kg
FoA = 170.34/(18.5*32+18.5*3.76*28);
B = 32.22*pi()/180;
pr1 = (((y+1)/2)*M*M*(sin(B))^2)
pr2 = (1+((y-1)/2)*M*M*(sin(B))^2)
pr3 = (((2*y)/(y+1))*M*M*(sin(B))^2-(y-1)/(y+1))
p_ratio_wedge = ((pr1/pr2)^(y/(y-1)))*(1/pr3)^(1/(y-1))
ploss_wedge = pa-pa*p_ratio_wedge
p_wedge = pa*p_ratio_wedge

pr4 = (1+((y-1)/2)*M*M)
pr5 = y*M*M*(sin(B))^2-((y-1)/2)
pr6 = M*M*(cos(B))^2
pr7 = pr2
M2 = sqrt(pr4/pr5+pr6/pr7)


pn1 = (((y+1)/2)*M2*M2);
pn2 = (1+((y-1)/2)*M2*M2);
pn3 = (((2*y)/(y+1))*M2*M2-(y-1)/(y+1));
p_ratio_nozzel = ((pn1/pn2)^(y/(y-1)))*(1/pn3)^(1/(y-1))
p_loss_nozzel = p_wedge-p_wedge*p_ratio_nozzel

total_loss = p_loss_nozzel+ploss_wedge % 8.2254e+03 Pa
total_ratio = (pa-total_loss)/pa       %0.5398

%The angle of the wedge can be increased, to induce a larger oblique shock 
%and therfore a smaller normal shock pressure loss. Two small shocks are
%better at maintaining stagnation pressure than a single big normal shock, 
%therefore the case with the bigger oblique shock has less total losses. 

%%%%%% b %%%%%%
%scramjet
%Atmosphere
p0a = pa*(1+(y-1)*M*M/2)^(y/(y-1));
T0a = Ta*(1+(y-1)*M*M/2);

%Diffuser/intake
p01 = p0a-total_loss;

tratio1 = ((1+((y-1)/2)*M*M)/(1+((y-1)/2)*M2*M2))
tratio2 = (1+((y-1)/2)*M*M*sin(B)^2)*(((2*y)/(y-1))*M*M*sin(B)^2)/((((y+1)^2)/(2*(y-1)))*M*M*sin(B)^2)
T01 = tratio1*tratio2*Ta
%Combustor
T04 = 1700
p04 = p01

phi = (T04-T01)*cp/(FoA*LHV)
%nozzel
Ve_c = sqrt(2*(y/(y-1))*R*T04.*(1-(pa./p04).^((y-1)/y)));

V = sos*M;
Thrust = (Ve_c-V) %(specific)  595
TSFC = 3600*4.44*phi*FoA./Thrust    % 0.7895 lbm/hr-lbf


%%%%%% c %%%%%%
%turbojet

HP_comp_r = 10;
while HP_comp_r<35;
%Diffuser
p01 = p0a;
T01 = T0a;


%HP compressor

p03 = HP_comp_r*p01;
n_comp = 0.9;
T03 = T01*(1+(1/n_comp)*(HP_comp_r^((y-1)/y)-1));

%Combustor
LHV = 43e6; %j/kg
% C12H26 + 18.5 (O2 + 3.76N2) -> 12 CO2 + 13 H2O + 69.56 N2
FoA = 170.34/(18.5*32+18.5*3.76*28);% fuel to air ratio
p04 = p03;


T04 = 1700;
phi = (T04-T03)*cp/(FoA*LHV);
%HP turbine
n_turb = 0.95;
T05 = T04 - (T01/n_comp)*(((HP_comp_r)^((y-1)/y))-1);
p05 = p04*(1-(1/(n_comp*n_turb*T04))*T01*(((HP_comp_r)^((y-1)/y))-1))^(y/(y-1));


%nozzel
Ve_c_t = sqrt(2*(y/(y-1))*R*T05.*(1-(pa./p05).^((y-1)/y)));


sos = sqrt(1.4*pa/rhoa);
V = sos*M;
Thrust = (Ve_c_t-V);
TSFC = 3600*4.44*phi*FoA./Thrust;
Np = (2*V/Ve_c_t)/(1+V/Ve_c_t);


figure(1)
subplot(1,3,1)
scatter(HP_comp_r,Thrust); hold on
xlabel('Compressor Ratio')
ylabel('Specific Thrust (N-s/kg)')
subplot(1,3,2)
scatter(HP_comp_r,Np); hold on
ylabel('Propulsive Efficiency')
xlabel('Compressor Ratio')
subplot(1,3,3)
scatter(HP_comp_r,TSFC); hold on
ylabel('TSFC (lb/hr-lb)')
xlabel('Compressor Ratio')

HP_comp_r = HP_comp_r + 3;

end

% Since the turbojet doesnt have a very big advantage of fuel consumption
% but is a lot worse at specific thrust, it means a turbojet would have to
% be much bigger than a scramjet in order to produce the same thrust. With
% this being said, size and drag is very important at supersonic speeds and
% so the scram jet would appear to be the smarter choice.
