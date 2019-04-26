###################################
%Material and System Variables
###################################
f = 8000/60; %Hz
omega = 2*pi*f; %Rad/s
engineV = 3000/100^3; %m^3
%TODO - Find and Justify choice of kinematic friction.
mu = .05;
%Ductile Iron Grade 65-45-12 
%http://www.matweb.com/search/DataSheet.aspx?MatGUID=0c178ad6ab2149d4a52a85a13645cc1a
ringMat = struct{'E',168*10^9,'v',.29,'Sy',332*10^6,'Sut',464*10^6,'rho',7.15*1000}; %Si Units
%310 Stainless Steel: High Creep Resistance and Temperature Resistance 
pistonMat = struct{'E',200*10^9,'v',.29,'Sy',205*10^6,'Sut',515*10^6,'rho',7750}; %Si Units
chamberMat = PistonMat;


########################
%Thermo Cycle Variables
########################

CR = 10;
%State 1
P_amb = 101325; %Pa
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
P3 = 8.187*10^9; %Pa
T3 = 2385;
%State 4
T4 = 1022;
u4 = 778;
rho4 = .777;
%Cycle Properties
w = (u3-u4)-(u2-u4)-p_amb(v_comp-v_amb);
m = engineV/(6*(v_4-v_comp));
Pow = 6*m*w*f/2;

#########
%Geometry
#########
%Piston
chamberDiameter = 7.88/100; %m
chamberArea = pi*chamberDiameter^2/4;
pistonSweep = (engineV/6)/(chamberArea);
pistonDiameter = chamberDiameter - 4/1000;
pistonArea = pi*pistonDiameter^2/4;
pistonHeadHeight = 50/1000;
pistonTailDiameter = 20/1000;
pistonTailArea = pi*pistonTailDiameter^2/4;
pistonTailHeight = 20/1000;
pistonVolume = pistonArea*pistonHeadHeight + pistonTailArea*pistonTailHeight;
pistonCoM = (pistonArea*pistonHeadHeight*(pistonHeadHeight/2) + pistonTailArea*pistonTailHeight*(pistonHeadHeight+pistonTailHeight/2))/pistonVolume;
%Rings
interference = 2/1000;
ringOuterDiameter = chamberDiameter + interference;
ringInnerDiameter = pistonDiameter - 10/1000;
ringHeight = 10/1000;
ringLocHeight = 15/1000;
ringLocHeight2 = pistonHeadHeight-15/1000;
%Connecting Rod
CRLength = 16/100; %m
%CrankShaft
crankRadius = pistonSweep/2; %m
%MaxTanAlpha
maxTan = crankRadius/sqrt(CRLength^2-crankRadius^2);
%TODO Ensure connecting rod doesn't collide with chamber.

##############################
%Piston and Ring Calculations
##############################
%Input Forces
maxPForce = P3*chamberArea;
maxReForce = maxPForce*maxTan;
%Spring Constant
k = chamberDiameter*ringHeight*ringMat.E/(.5*(chamberDiameter-ringInnerDiameter));
%Solving Geometric Constants
alpha = pistonHeadHeight*(2*pistonCoM-ringLocHeight-ringLocHeight2) + pistonTailHeight*(2*pistonCoM-ringLocHeight-ringLocHeight2);
beta = -4*pistonCoM^2+3*pistonCoM*(ringLocHeight+ringLocHeight2)-ringLocHeight^2-ringLocHeight2^2;
denm = (2*k*(ringLocHeight^2-2*ringLocHeight*ringLocHeight2+ringLocHeight2^2));
x = maxReForce*(alpha+beta)/denm;
theta = maxReForce*(2*pistonHeadHeight+2*pistonTailHeight-4*pistonCoM + ringLocHeight + ringLocHeight2)/denm;
%Ring Forces
ringForces = [k*(x-(pistonCoM-ringLocHeight)*theta),k*(-x+(pistonCoM-ringLocHeight)*theta),k*(x-(pistonCoM-ringLocHeight2)*theta),k*(-x+(pistonCoM-ringLocHeight2)*theta)];
%Ring Pressures
minPFit = P3 - min(ringForces)/(chamberDiameter*ringHeight);

denm = (chamberDiameter/pistonMat.E)*(1+pistonMat.v) + (ringOuterDiameter/ringMat.E)*(1+ringMat.v);
PFit = interference/denm;
%Factor Of Safety
if 1.5*minPFit > PFit
   error('Ring Pressure does not meet mininum condition.')
endif

%Stress
maxRingRadialStress = PFit + max(ringForces)/(chamberDiameter*ringHeight);
maxRingTangentStress = maxRingRadialStress;
maxRingShearStress = 3*(chamberDiameter^2-pistonDiameter^2)*P3/(4*ringHeight*pistonDiameter) + 3*mu*PFit*chamberDiameter^2/pistonDiameter^2;
ringVonMises = sqrt(maxRingRadialStress^2 + 3*maxRingShearStress^2);
if 3*ringVonMises > ringMat.Sy
  error("Ring stresses exceed factor of safety.")
endif
