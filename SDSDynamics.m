% SDS Propagation Method
function dSSTTpropagation = SDSDynamics(t,S,constants)
l = wrapToPi(S(1)); g = wrapToPi(S(2)); h = wrapToPi(S(3)); L = S(4); G = S(5); H = S(6); ks = S(7); Ks = S(8);

SunCOE = getSun(constants,t*3600);
SunDel = normalize(COE_to_Delaunay(SunCOE,constants.MU,1),constants.rE,'vec','Del',1); % Using normalized mu for sun propagation, we only care about h from the Sun coordinates.
hs = SunDel(3); % rad
cs = SunDel(6) / SunDel(5);
sSun = sqrt(1 - cs^2);
nu = constants.MU_units^2 / SunDel(4)^3;

% Using normalized constants for propagation
RE = 1;
J2 = constants.j2;
mu = constants.mu_units;
Aom = constants.Aom / constants.rE^2; % km^2/kg 
Psrp = constants.Psrp * constants.rE * 3600^2; % kg/s^2/km
rho = constants.rho;
b = (1+rho) * Aom * Psrp;

dS = MeanDynamicsFunction(G,H,J2,L,RE,b,cs,g,h,hs,ks,mu,nu,sSun);

dim = constants.STTDim;
STTdot1 = getJOrder1(G,H,J2,L,RE,b,cs,g,h,hs,ks,mu,sSun);
STM1 = reshape(S(dim+1:(dim+dim^2)),dim,dim)';
phiDot = STTdot1 * STM1;
phiDotVec = reshape(phiDot',dim^2,1);

if constants.STTOrder == 2
    STTdot2 = getJOrder2(G,H,J2,L,RE,b,cs,g,h,hs,ks,mu,sSun);
    phiDotVec = [phiDotVec; Compute2ndOrderSTTdot(S((dim+dim^2+1):end), STTdot2, dim, STTdot1, STM1)];
end

dSSTTpropagation = [dS;phiDotVec];

end
