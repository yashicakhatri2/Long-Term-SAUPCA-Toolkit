function dSFull2 = EOM(t, S, constants)
rvec = S(1:3); r = norm(rvec);
rI = rvec(1); rJ = rvec(2); rK = rvec(3);
mu = constants.mu;
rE = constants.rE;
J2 = constants.j2;

% 2BP
a2BP = - mu / r^3 * rvec;

% SRP
SSEvec = COE_to_Cartesian_2(getSun(constants,t),constants.MU,1);
rrelvec = SSEvec(1:3) - rvec; rrel = norm(rrelvec);
b = (1+constants.rho) * constants.Psrp * constants.Aom;
aSRP = -b * rrelvec / rrel;

% J2
aI = -3 * J2 * mu * rE^2 * rI / 2 / r^5 * (1 - 5 * rK^2 / r^2);
aJ = -3 * J2 * mu * rE^2 * rJ / 2 / r^5 * (1 - 5 * rK^2 / r^2);
aK = -3 * J2 * mu * rE^2 * rK / 2 / r^5 * (3 - 5 * rK^2 / r^2);
aJ2 = [aI;aJ;aK];
% aJ2 = ECEFtoECI(aJ2_,t,1);

atot = a2BP + aSRP + aJ2;

dSFull2(1:3,1) = S(4:6);
dSFull2(4:6,1) = atot;

end