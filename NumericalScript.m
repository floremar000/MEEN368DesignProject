
function NumericalScript
  ###################################
  %Material and System Variables
  ###################################
  f = 8000/60; %Hz
  global omega = 2*pi*f; %Rad/s
  engineV = 3000/100^3; %m^3
  period = 2/f;
  gamma = 1.4;
  %TODO - Find and Justify choice of kinematic friction.
  mu = .05;

  %Ductile Iron Grade 65-45-12 
  %http://www.matweb.com/search/DataSheet.aspx?MatGUID=0c178ad6ab2149d4a52a85a13645cc1a
  %https://www.ductile.org/didata/Section3/3part1.htm#Fatigue%20Limit
  ringMat = struct('E',168*10^9,'v',.29,'Sy',332*10^6,'Sut',464*10^6,'SeN',117*10^6,'rho',7.15*1000); %Si Units
  %Titanium Ti-6Al-4V 
  pistonMat = struct('E',113.8*10^9,'v',.29,'Sy',880*10^6,'Sut',950*10^6,'Se',510*10^6,'rho',4430); %Si Units
  chamberMat = pistonMat;
  crankShaftMat = struct('E',200*10^9,'v',.29,'Sy',580*10^6,'Sut',690*10^6,'rho',7850); %Si Units

  ########################
  %Thermo Cycle Variables
  ########################

  CR = 10;
  %State 1
  P_amb = 101325; %kPa
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
  P3 = 8.187*10^6; %Pa
  T3 = 2385;
  %State 4
  T4 = 1022;
  u4 = 778;
  rho4 = .777;
  v_4 = 1/rho4;
  %Cycle Properties
  w = (u3-u4)-(u2-u_amb)-P_amb*(v_comp-v_amb)/1000;
  airMass = engineV/(6*(v_4-v_comp));
  Pow = 6*airMass*w*f/2;
  
  %Volumes
  V_comp = v_comp*airMass;
  V_4 = v_4*airMass;
  
  #########
  %Geometry
  #########
  %Piston
  chamberDiameter = 7.88/100; %m
  chamberArea = pi*chamberDiameter^2/4;
  pistonSweep = (engineV/6)/(chamberArea);
  %Piston Head
  pistonHeadDiameter = chamberDiameter - 4/1000;
  pistonHeadArea = pi*pistonHeadDiameter^2/4;
  pistonHeadHeight = 65/1000;
  %Piston Tail
  pistonTailDiameter = 40/1000;
  pistonTailArea = pi*pistonTailDiameter^2/4;
  pistonTailHeight = 20/1000;
  pistonTailFillet = 8/1000;
  pistonTailStressConcentration = 1.5;
  %Joint Hole
  pistonJointDiameter = 15/1000;
  if pistonJointDiameter>=pistonTailHeight
    error('Piston Joint Hole does not fit tail.')
  endif
  pistonForceHeight = pistonTailHeight/2;
  angle = 2*acos(pistonJointDiameter/pistonTailDiameter);
  pistonChordArea = pistonTailDiameter^2*(angle - sin(angle))/8;
  pistonVolume = pistonHeadArea*pistonHeadHeight - .8^2*pistonHeadArea*pistonHeadHeight + pistonTailArea*pistonTailHeight;
  pistonCoM = (pistonHeadArea*pistonHeadHeight*(pistonHeadHeight/2) + pistonTailArea*pistonTailHeight*(pistonHeadHeight+pistonTailHeight/2))/pistonVolume;
  pistonMass = pistonMat.rho*pistonVolume;
  %Piston Channel
  pistonChannelDepth = 5/1000;
  pistonRingChannelFillet = 2/1000;
  pistonChannelDiameter = pistonHeadDiameter - 2*pistonChannelDepth;
  %Stress Concentration in Channels
  %Round shaft with flat-bottom groove in bending and/or tension.
  pistonChannelConcentrationFactor = 3.5;
  numPistons = 6;
  %Intake Hole
  intakeArea = chamberArea/5;

  %Rings
  interference = .0009*chamberDiameter;
  ringOuterDiameter = 7.89/100;
  ringInnerDiameter = pistonChannelDiameter;
  ringHeight = 5/1000;
  ringLocHeight = 15/1000;
  ringLocHeight2 = pistonHeadHeight-15/1000;

  %Connecting Rod
  global CRLength = 16/100; %m
  CRMass = 2;
  
  %CrankShaft Geometry
  global crankRadius = pistonSweep/2; %m
  crankLength = pistonHeadDiameter;
  shaftLength = 10/1000;
  crankConnectorThickness = 20/1000;
  crankShaftSubsectionLength = crankLength + shaftLength + 2*crankConnectorThickness;
  shaftBearingDistance = 30/1000;
  crankShaftLength = shaftBearingDistance +crankShaftSubsectionLength*6;
  shaftDiameter = 75/1000;
  crankDiameter = shaftDiameter;
  shaftI = pi*shaftDiameter^4/(4*2^4);
  crankI  = pi*crankDiameter^4/(4*2^4);
  crankConnectorI = shaftDiameter*(crankRadius+shaftDiameter/2+crankDiameter/2)^3/12;
  %MaxTanAlpha
  maxTan = crankRadius/sqrt(CRLength^2-crankRadius^2);
  %TODO Ensure connecting rod doesn't collide with chamber.
  
  
  ##############################
  %Forces
  ##############################  
  samples = 10;
  t = 0:period/samples:period;
  global initialPistonHeight = V_comp/chamberArea;
  
  tguess = 1.5*pi/omega;
  criticalTime = fzero(@pistonPositionSol,tguess);
  criticalTime = mod(criticalTime,2*pi);
  if omega*criticalTime <= pi
    error('Wrong Time Solution')
  endif
  timeSamples = length(t);
  pistonPosition= zeros(numPistons,timeSamples);
  pistonVelocity= zeros(numPistons,timeSamples);
  pistonAcceleration= zeros(numPistons,timeSamples);
  gasForceMagnitude= zeros(numPistons,timeSamples);
  totalForceMagnitude =zeros(numPistons,timeSamples);
  tn = t;
  for i = 2:numPistons
    t = [t;t(1,:)+period*(i-1)/(numPistons)];
  endfor

  %Piston Position
  pistonPosition = crankRadius + CRLength - (crankRadius.*cos(omega.*t) + sqrt(CRLength^2-(crankRadius*sin(omega*t)).^2));
  %Piston Velocty
  denm = sqrt(CRLength^2-crankRadius^2.*sin(omega.*t).^2);
  pistonVelocity = omega*crankRadius^2.*sin(omega*t).*cos(omega*t)./denm + omega*crankRadius*sin(omega*t);
  %Piston Acceleration
  denm1 = sqrt(CRLength^2-crankRadius^2.*sin(omega.*t).^2);
  denm2 = (CRLength^2-crankRadius^2.*sin(omega.*t).^2).^(1.5);
  pistonAcceleration = omega^2*crankRadius^4*(sin(omega*t).*cos(omega*t)).^2./denm2 + omega^2*crankRadius^2*(cos(omega*t).^2-sin(omega*t).^2)./denm1 + omega^2*crankRadius*cos(omega*t);
  
  
  disp('Starting Force Calculations')
  %Stroke Timing
  powerStroke = (mod(omega*t,4*pi) <= pi)&(mod(omega*t,4*pi) >= 0);
  exhaustStroke = ((mod(omega*t,4*pi) <= 2*pi)&(mod(omega*t,4*pi) > pi))|((mod(omega*t,4*pi) <= 2*pi + omega*criticalTime)&(mod(omega*t,4*pi) > 3*pi));
  intakeStroke = ((mod(omega*t,4*pi) <= 3*pi)&(mod(omega*t,4*pi) > 2*pi));
  compressionStroke = ((mod(omega*t,4*pi) < 4*pi)&(mod(omega*t,4*pi) > 2*pi+omega*criticalTime));
  %Gas Force Calculations
  fPower = chamberArea*(P3*(V_comp./(chamberArea*(pistonPosition+initialPistonHeight))).^gamma - P_amb);
  fPower = fPower.*powerStroke;
    
  fExhaust = .5*rho_amb*chamberArea.*(chamberArea*pistonVelocity/intakeArea).^2;
  fExhaust = fExhaust.*exhaustStroke;
    
  fIntake = -.5*rho_amb*chamberArea.*(chamberArea*pistonVelocity/intakeArea).^2;
  fIntake = fIntake.*intakeStroke;
    
  fCompression = chamberArea.*(P_amb*(10*initialPistonHeight./(pistonPosition+initialPistonHeight)).^gamma - P_amb);
  fCompression = fCompression.*compressionStroke;
  gasForceMagnitude = fPower+fExhaust+fIntake+fCompression;
  totalForceMagnitude = gasForceMagnitude-(pistonMass+CRMass)*pistonAcceleration;
  totalPistonForceMagnitude = gasForceMagnitude-(pistonMass)*pistonAcceleration;
  %Rows are time vector
  %Columns are different pistons
  %3rd is vecotrs i and j.
  tanAlpha = crankRadius*sin(omega*t)./sqrt(CRLength^2-(crankRadius*sin(omega*t)).^2);
  totalForceVectors = zeros([size(t),3]);
  totalForceVectors(:,:,1) = -totalForceMagnitude;
  totalPistonForceVectors(:,:,1) = -totalPistonForceMagnitude;
  totalForceVectors(:,:,2) = tanAlpha.*totalForceMagnitude;
  totalPistonForceVectors(:,:,2) = tanAlpha.*totalPistonForceMagnitude;
  
  plot(t(1,:),totalForceVectors(:,:,1));
  
  
  ##############################
  %Ring Calculations
  ##############################
  %Input Forces
  maxPForce = max(max(abs(totalPistonForceVectors(:,:,1))));
  maxReForce = max(max(abs(totalPistonForceVectors(:,:,2))));
  %Spring Constant
  k = chamberDiameter*ringHeight*ringMat.E/(.5*(chamberDiameter-ringInnerDiameter));
  %Solving Geometric Constants
  alpha = pistonHeadHeight*(2*pistonCoM-ringLocHeight-ringLocHeight2) + pistonForceHeight*(2*pistonCoM-ringLocHeight-ringLocHeight2);
  beta = -4*pistonCoM^2+3*pistonCoM*(ringLocHeight+ringLocHeight2)-ringLocHeight^2-ringLocHeight2^2;
  denm = (2*k*(ringLocHeight^2-2*ringLocHeight*ringLocHeight2+ringLocHeight2^2));
  x = maxReForce*(alpha+beta)/denm;
  theta = maxReForce*(2*pistonHeadHeight+2*pistonForceHeight-4*pistonCoM + ringLocHeight + ringLocHeight2)/denm;
  %Ring Forces
  ringForces = [k*(x-(pistonCoM-ringLocHeight)*theta),k*(-x+(pistonCoM-ringLocHeight)*theta),k*(x-(pistonCoM-ringLocHeight2)*theta),k*(-x+(pistonCoM-ringLocHeight2)*theta)];
  %Ring Pressures
  minPFit = P3 - min(ringForces)/(chamberDiameter*ringHeight);

  denm = (chamberDiameter/(200*10^9))*(1+pistonMat.v) + (ringOuterDiameter/ringMat.E)*(1+ringMat.v);
  PFit = interference/denm;

  %Factor Of Safety
  if 1.5*minPFit > PFit
     error('Ring Pressure does not meet mininum condition.')
  endif

  %Ring Stresses
  maxRingRadialStress = PFit + max(ringForces)/(chamberDiameter*ringHeight);
  maxRingTangentStress = maxRingRadialStress;
  maxRingShearStress = 3*(chamberDiameter^2-pistonHeadDiameter^2)*P3/(4*ringHeight*pistonHeadDiameter) + 3*mu*PFit*chamberDiameter^2/pistonHeadDiameter^2;
  maxRingVonMises = sqrt(maxRingRadialStress^2 + 3*maxRingShearStress^2);

  minRingRadialStress = PFit + min(ringForces)/(chamberDiameter*ringHeight);
  minRingTangentStress = minRingRadialStress;
  minRingShearStress = -3*mu*PFit*chamberDiameter^2/pistonHeadDiameter^2;
  ampRingRadialStress = (maxRingRadialStress-minRingRadialStress)/2;
  ampRingRadialStress = (maxRingRadialStress-minRingRadialStress)/2;
  ampRingTangentStress = ampRingRadialStress;
  %Stresses for fatiude Criteria
  ampRingVonMises = sqrt((ampRingRadialStress/.85)^2 +3*ampRingTangentStress^2);
  %Since the mean stress radial stress is compressive, it's igorned.
  meanRingShearStress = (maxRingShearStress+minRingShearStress)/2;
  meanRingVonMises = sqrt(3*meanRingShearStress^2);
  %mod-Goodman modified Formula
  nRingFatigue = inv(meanRingVonMises/ringMat.Sy + ampRingVonMises/ringMat.SeN);

  if nRingFatigue < 1.5
    disp('Calulated n');
    disp(nRingFatigue);
    error('Ring stresses do not meet Fatigue Criteria')
  endif


  if 3*maxRingVonMises > ringMat.Sy
    disp('Ratio of VonMises to Yield')
    disp(maxRingVonMises/ringMat.Sy)
    error("Ring stresses exceed factor of safety.")
  endif
  
  
  ###############################
  %Piston Stresses and Deflection
  ###############################
  pistonHeadStress = -maxPForce/pistonHeadArea;
  pistonChannelPressureStress = pistonHeadStress*pistonHeadDiameter^2/(pistonChannelDiameter^2);
  pistonChannelPushoutStress = -P3*(ringOuterDiameter^2-pistonHeadDiameter^2)/(pistonHeadDiameter^2-ringInnerDiameter^2);
  pistonChannelMaxStress = pistonChannelConcentrationFactor*(pistonChannelPressureStress+pistonChannelPushoutStress);

  pistonTailStress = -maxPForce/pistonTailArea;
  pistonTailMaxStress = pistonTailStressConcentration*pistonTailStress;

  pistonJointShearStress = 3*maxReForce/(4*pistonChordArea);
  pistonJointNormalStress = -maxPForce/(pistonTailArea-2*pistonChordArea);
  pistonJointPulloutStress = maxPForce/(2*pistonChordArea);

  [pistonMaxVonMises, maxindex] = max([abs(pistonJointNormalStress),abs(sqrt(3)*pistonJointShearStress),abs(pistonTailMaxStress), abs(pistonChannelMaxStress),abs(pistonJointPulloutStress)]);
  %Axial go from roughly 0 to max. Mean stress is compressive.
  pistonAlternatingAxialVonMises = max([abs(pistonJointNormalStress)/.85,abs(pistonTailMaxStress)/.85, abs(pistonChannelMaxStress)/.85]);
  %Push Out Fatigue
  pistonMeanPushOutVonMises = sqrt(3)*pistonJointShearStress/2;
  pistonAmpPushOutVonMises = sqrt(3)*pistonJointShearStress/2;

  if 3*pistonMaxVonMises > pistonMat.Sy
    disp('Ratio of VonMises to Yield')
    disp(pistonMaxVonMises/pistonMat.Sy)
    disp(maxindex)
    error("Piston stresses exceed factor of safety.")
  endif
  %Axial Fatigue
  if 2*pistonAlternatingAxialVonMises > pistonMat.Sut/2
    disp('Ratio of Axial Von Mises to Se')
    disp(2*pistonAlternatingAxialVonMises/pistonMat.Sut)
    error('Alternating axial piston stresses exceed factor of safety')
  endif
  %Push out Fatigue
  pistonPushOutFatigueFoS = inv(pistonMeanPushOutVonMises/pistonMat.Sut + 2*pistonAmpPushOutVonMises/pistonMat.Sut);
  if pistonPushOutFatigueFoS < 1.5
    disp('Ratio of Push Out Von Mises to Se')
    disp(2*pistonAmpPushOutVonMises/pistonMat.Sut)
    error('Alternating piston shear stresses exceed factor of safety')
  endif

  %Estimated Piston Deflection
  pistonHeadDeflection = pistonHeadHeight*pistonHeadStress/pistonMat.E;
  pistonTailDeflection = pistonForceHeight*pistonTailStress/pistonMat.E;
  pistonTotalDeflection = pistonHeadDeflection + pistonTailDeflection;

  
  ###############################
  %CrankShaft Calculations
  ###############################
  radiusVectors = zeros([size(t),3]);
  radiusVectors(:,:,1) = crankRadius*cos(omega*t);
  radiusVectors(:,:,2) = crankRadius*sin(omega*t);
  

  R1 = zeros([1,timeSamples,3]);
  R2 = zeros([1,timeSamples,3]);
  
  for i = 1:numPistons
    R2 = R2-totalForceVectors(i,:,:)*(crankShaftSubsectionLength*(i-1)+shaftBearingDistance+crankLength/2+crankConnectorThickness)/crankShaftLength;
  endfor
  R1 = -R2;
  for i = 1:numPistons
    R1 = R1 - totalForceVectors(i,:,:);
  endfor
  
  crankShaftLocations = 0:crankShaftLength/samples:crankShaftLength;
  numLocations = length(crankShaftLocations);
  %Crank Shaft Internal Reaction Equations
  crankShaftTorque = zeros([1,timeSamples,3,numLocations]);
  for i = 1:numLocations
    loc = crankShaftLocations(i);
    for j = 1:numPistons
      crankShaftTorque(1,:,:,i) = crankShaftTorque(1,:,:,i) + cross(totalForceVectors(j,:,:),radiusVectors(j,:,:),3)*heaviside(loc-crankShaftSubsectionLength*(j-1)-shaftBearingDistance-crankLength/2-crankConnectorThickness);
    endfor
  endfor
  
  
  crankShaftShear = zeros([1,timeSamples,3,numLocations]);
  for i = 1:numLocations
    loc = crankShaftLocations(i);
    crankShaftShear(1,:,:,i) = R1*heaviside(loc);
    for j = 1:numPistons
      crankShaftShear(1,:,:,i) = crankShaftShear(1,:,:,i) + totalForceVectors(j,:,:)*heaviside(loc-crankShaftSubsectionLength*(j-1)-shaftBearingDistance-crankLength/2-crankConnectorThickness);
    endfor
    crankShaftShear(1,:,:,i) = crankShaftShear(1,:,:,i)  + R2*heaviside(loc-crankShaftLength);
  endfor
  
  crankShaftMoment = zeros([1,timeSamples,3,numLocations]);
  for i = 1:numLocations
    loc = crankShaftLocations(i);
    crankShaftMoment(1,:,:,i) = R1*heaviside(loc)*loc;
    for j = 1:numPistons
      arm = loc-crankShaftSubsectionLength*(j-1)-shaftBearingDistance-crankLength/2-crankConnectorThickness;
      crankShaftMoment(1,:,:,i) = crankShaftMoment(1,:,:,i) + totalForceVectors(j,:,:)*heaviside(arm)*arm;
    endfor
    crankShaftMoment(1,:,:,i) = crankShaftMoment(1,:,:,i)  + R2*heaviside(loc-crankShaftLength)*(loc-crankShaftLength);
  endfor
  
  crankShaftSlope = zeros([1,timeSamples,3,numLocations]);
  for i = 1:numLocations
    loc = crankShaftLocations(i);
    crankShaftSlope(1,:,:,i) = .5*R1*heaviside(loc)*loc.^2;
    for j = 1:numPistons
      arm = loc-crankShaftSubsectionLength*(j-1)-shaftBearingDistance-crankLength/2-crankConnectorThickness;
      crankShaftSlope(1,:,:,i) = crankShaftSlope(1,:,:,i) + .5*totalForceVectors(j,:,:)*heaviside(arm)*arm.^2;
    endfor
  endfor

    
  crankShaftDisplacement = zeros([1,timeSamples,3,numLocations]);
  for i = 1:numLocations
    loc = crankShaftLocations(i);
    crankShaftDisplacement(1,:,:,i) = R1*heaviside(loc)*loc.^3/6;
    for j = 1:numPistons
      arm = loc-crankShaftSubsectionLength*(j-1)-shaftBearingDistance-crankLength/2-crankConnectorThickness;
      crankShaftDisplacement(1,:,:,i) = crankShaftDisplacement(1,:,:,i) + totalForceVectors(j,:,:)*heaviside(arm)*arm.^3/6;
    endfor
  endfor
  
  C1 = -crankShaftDisplacement(1,:,:,end)/crankShaftLength;
  for i = 1:numLocations
    loc = crankShaftLocations(i);
    crankShaftSlope(1,:,:,i) = crankShaftSlope(1,:,:,i) + C1;
    crankShaftDisplacement(1,:,:,i) = crankShaftDisplacement(1,:,:,i) + C1*loc;
  endfor
  crankShaftSlope = crankShaftSlope/(pistonMat.E*shaftI);
  crankShaftDisplacement = crankShaftDisplacement/(pistonMat.E*shaftI);
  
  crankShaftDisplacementMag = crankShaftDisplacement(1,:,1,:).^2+crankShaftDisplacement(1,:,2,:).^2;
  crankShaftDisplacementMag = sqrt(crankShaftDisplacementMag);
  crankShaftDisplacementMax = max(max(crankShaftDisplacementMag));
  
  crankShaftSlopeMag = crankShaftSlope(1,:,1,:).^2+crankShaftSlope(1,:,2,:).^2;
  crankShaftSlopeMag = sqrt(crankShaftSlopeMag);
  crankShaftSlopeMax = max(max(crankShaftSlopeMag));
  
  crankShaftMomentMag = crankShaftMoment(1,:,1,:).^2 + crankShaftMoment(1,:,2,:).^2;
  crankShaftMomentMag = sqrt(crankShaftMomentMag);
  maxShaftMoment = max(max(abs(crankShaftMomentMag)));
  crankShaftMomentStress = maxShaftMoment*shaftDiameter/(2*shaftI);
  
  crankShaftTorqueMax = max(max(max(abs(crankShaftTorque))));
  crankShaftShearStress = crankShaftTorqueMax*shaftDiameter/(4*shaftI);
  
  crankShaftVonMises = sqrt(crankShaftMomentStress^2+3*crankShaftShearStress^2);
  
  if 1.5*crankShaftVonMises > crankShaftMat.Sy
    error('Crank Shaft fails to meet FoS')
  endif
  
  
  
  
  
  disp('Paused')
 save vars.mat
endfunction

function pisPos = pistonPositionSol(t)
  global crankRadius CRLength omega initialPistonHeight
  pisPos = crankRadius + CRLength - (crankRadius.*cos(omega.*t) + sqrt(CRLength^2-crankRadius^2.*sin(omega.*t).^2)) - 9*initialPistonHeight;
endfunction
 
 