clc; clear all; close all;
M = 200/761;
pa = 101350% Pa
rhoa = 1.225% kg/m3
Ta = 288% K
V = 200*0.44704
y = 1.4
R = 287; %j/kg-K
cp = y*R/(y-1);

%Setting fixed parameters
HP_comp_r = 40
T04 = linspace(1500,2500,1000); 
BPR = 13;

p0a = pa*(1+(y-1)*M*M/2)^(y/(y-1));
T0a = Ta*(1+(y-1)*M*M/2);
p01 = p0a;
T01 = T0a;
fan_r = 1.5;
p02 = fan_r*p01;
n_fan = 0.9;
T02 = T01*(1+(1/n_fan)*(fan_r^((y-1)/y)-1));
p03 = HP_comp_r*p02;
n_comp = 0.9;
T03 = T02*(1+(1/n_comp)*(HP_comp_r^((y-1)/y)-1));
LHV = 43e6; %j/kg
FoA = 170.34/(18.5*32+18.5*3.76*28);% fuel to air ratio
p04 = p03;
phi = (T04-T03)*cp./(FoA*LHV)
n_turb = 0.95;
T05 = T04 - (T02./n_comp).*(((HP_comp_r).^((y-1)/y))-1);
p05 = p04.*(1-(1./(n_comp.*n_turb.*T04)).*T02.*(((HP_comp_r).^((y-1)/y))-1)).^(y/(y-1));
T06 = T05 - (1+BPR).*(T01./n_fan).*(((fan_r).^((y-1)/y))-1);
p06 = p05.*(1-((1+BPR)./(n_fan.*n_turb.*T05)).*T01.*(((fan_r).^((y-1)/y))-1)).^(y/(y-1));
Ve_c = sqrt(2.*(y/(y-1))*R*T06.*(1-(pa./p06).^((y-1)/y)));
Ve_f = sqrt(2.*(y/(y-1))*R*T02.*(1-(pa./p02).^((y-1)/y)));
Thrust = (BPR*(Ve_f)+(Ve_c)) %Gross thurst used instead
Np = Thrust.*V./(BPR*((Ve_f.^2)/2-(V.^2)/2)+((Ve_c.^2)/2-(V.^2)/2))
TSFC = 3600*4.44*2.2*phi*FoA./Thrust


figure(1)
subplot(1,2,1)
plot(T04,Thrust); hold on
xlabel('to4')
ylabel('Specific Thrust (N-s/kg)')
subplot(1,2,2)
plot(T04,phi); hold on
ylabel('Phi')
xlabel('T04')
