###################################
%Material and System Variables
###################################
f = 8000/60; %Hz
omega = 2*pi*f; %Rad/s
EngineV = 3000/100^3; %m^3

%Ductile Iron Grade 65-40-18
GasketMat = struct{'E',205*10^9,'v',.29,'Sy',370*10^6,'Sut',440*10^6,'rho',7.87*1000}; %Si Units
%310 Stainless Steel: High Creep Resistance and Temperature Resistance 
PistonMat = struct{'E',200*10^9,'v',.29,'Sy',205*10^6,'Sut',515*10^6,'rho',7750}; %Si Units
ChamberMat = PistonMat;


########################
%Thermo Cycle Variables
########################

CR = 10;
%State 1
p_amb = 101325; %Pa
T_amb = 288; %K
rho_amb = 1.225; %kg/m^3
v_amb = 1/rho_amb; %m^3/kg
u_amb = 205.71; %kJ/kg
%State 2
v_comp = v_amb/10;
rho_comp = 10*rho_amb;
T2 = 703; %K
u2 = 514.95;
%State 3
u3 = 2051.1; %kJ/kg
p3 = 8.187*10^9; %Pa
T3 = 2385;
%State 4
T4 = 1022;
u4 = 778;
rho4 = .777;
%Cycle Properties
w = (u3-u4)-(u2-u4)-p_amb(v_comp-v_amb);
m = EngineV/(6*(v_4-v_comp));
P = 6*m*w*f/2;

#########
%Geometry
#########
%Piston
ChamberDiameter = 7.88/100; %m
ChamberArea = pi*ChamberDiameter^2/4;
PistonSweep = (EngineV/6)/(ChamberArea);
%Gaskets
%Connecting Rod
CRLength = 16/100; %m