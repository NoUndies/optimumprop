function [C BETA] = optimumprop(aerofunc, BLADE, cd0, vel, rpm, power, blades, rho)

  % This function uses Goldstein's equations in conjunction with BEMT to produce 
  % an ideal prop geometry. See Design of Optimum Propellers, Adkins & Liebeck, 1994.
  % -- Sean Morrison, 2017

  element = BLADE(end)-BLADE(end-1);
  rad = BLADE(end)+element/2;
  omega = rpm*2*pi/60
  zeta = 0;
  Cp = 2*power/rho/vel^3/pi/rad^2;
  XI = BLADE/rad;
  lambda = vel/omega/rad;
  X = omega*BLADE/vel;
  a = 0.5;
  err = 1;

  ALF = 0:1:90;
  ALF = ALF*pi/180;
  [CL_ CD_] = aerofunc(ALF);
  LND = CL_./(CD_+cd0);
  m = max(LND(:));
  idx = find(LND==m);
  ALPHA = ALF(idx)*ones(size(BLADE));
  CL = CL_(idx);
  EPSILON = 1/LND(idx);

  while err > 1e-5
    phi_t = atan(lambda*(1+zeta/2));
    PHI = atan(tan(phi_t)./XI);
    f = blades/2*(1-XI)/sin(phi_t);
    F = 2/pi*acos(exp(-f));
    G = F.*cos(PHI).*sin(PHI);
    WC = 4*pi*lambda.*G*vel*rad*zeta/blades./CL;
    A = zeta/2./X.*(cos(PHI)).^2.*(1-EPSILON.*tan(PHI));
    B = zeta/2*(cos(PHI)).*sin(PHI).*(1+EPSILON./tan(PHI));
    W = vel*(1+A)./sin(PHI);
    C = WC./W;
    BETA = ALPHA+PHI;
    I1p = 4*XI.*G.*(1-EPSILON.*tan(PHI));
    I2p = lambda.*I1p/2./XI.*(1+EPSILON./tan(PHI)).*sin(PHI).*cos(PHI);
    J1p = 4*XI.*G.*(1+EPSILON./tan(PHI));
    J2p = (J1p/2).*(1-EPSILON.*tan(PHI)).*cos(PHI).*cos(PHI); 
    I1 = sum(I1p);
    I2 = sum(I2p);
    J1 = sum(J1p);
    J2 = sum(J2p);
    tmp = J1/2./J2-sqrt((J1/2./J2).^2-Cp./J2);
    err = abs(zeta-tmp)
    zeta = (1-a)*zeta+tmp*a
  end

  Y0s = 0.*cos(BETA);
  Y1s = -C.*cos(BETA);
  Z0s = 0.*sin(BETA);
  Z1s = -C.*sin(BETA);
  Xs = BLADE;

  for i=1:size(BLADE,2)-1
    x = [Xs(i) Xs(i)];
    y = [Y0s(i) Y1s(i)];
    z = [Z0s(i) Z1s(i)];
    plot3(x,y,z);
    hold on
  end

  plot3(Xs,Y0s,Z0s)
  hold on
  plot3(Xs,Y1s,Z1s)
  axis([0 1 -0.5 0.5 -0.5 0.5]);
end
