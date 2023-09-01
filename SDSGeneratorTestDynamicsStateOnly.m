% SDS Propagation Method
function dS = SDSGeneratorTestDynamicsStateOnly(t,S,constants)
l = wrapToPi(S(1)); g = wrapToPi(S(2)); h = wrapToPi(S(3)); L = S(4); G = S(5); H = S(6); ks = wrapToPi(S(7)); Ks = S(8);

SunCOE = getSun(constants,t*3600);
SunCart = COE_to_Cartesian_2(SunCOE, constants.MU,1);
SunDel = normalize(COE_to_Delaunay(SunCOE,constants.MU,1),constants.rE,'vec','Del',1); % Using normalized mu for sun propagation, we only care about h from the Sun coordinates.
SatCOE = COE_to_Delaunay(normalize([l;g;h;L;G;H],constants.rE,'vec','Del',0),constants.mu,0);
SatCart = COE_to_Cartesian_2(SatCOE, constants.mu,1);
hs = SunDel(3); % rad
cs = SunDel(6) / SunDel(5);
sSun = sqrt(1 - cs^2);
nu = constants.MU_units^2 / SunDel(4)^3;

% Using normalized constants for propagation
E = 1;
RE = 1;
J2 = constants.j2;
mu = constants.mu_units;
Aom = constants.Aom / constants.rE^2; % km^2/kg 
Psrp = constants.Psrp * constants.rE * 3600^2; % kg/s^2/km
rho = constants.rho;

b = (1+rho) * Aom * Psrp;

a = L^2/mu;
e = sqrt(1-(G/L)^2);
[Ecc, TA] = MAtoEccTA(l,e);
rss = norm(SunCart(1:3) - SatCart(1:3)) / constants.rE;
r = a * (1 - e * cos(Ecc));

dS = MeanDynamicsFunction(G,H,J2,L,RE,b,cs,g,h,hs,ks,mu,nu,sSun) + ShortPeriodDynamicsFunction(Ecc,G,H,J2,L,RE,TA,b,cs,g,h,hs,ks,mu,r,sSun); 
% dS = HtildeDynamicsFunction(G,H,J2,L,RE,TA,b,cs,g,h,hs,ks,mu,nu,r,rss,sSun);

end

