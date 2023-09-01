% SDS Propagation Method
function dSSTTpropagation = SDSDynamicsForLongTCA(t,S,constants)
l1 = wrapToPi(S(1)); g1 = wrapToPi(S(2)); h1 = wrapToPi(S(3)); L1 = S(4); G1 = S(5); H1 = S(6); ks1 = S(7); Ks1 = S(8);
l2 = wrapToPi(S(585)); g2 = wrapToPi(S(586)); h2 = wrapToPi(S(587)); L2 = S(588); G2 = S(589); H2 = S(590); ks2 = S(591); Ks2 = S(592);
% l2 = l1 -dl; g2 = g1 - dg; h2 = h1 - dh; L2 = L1 - dL; G2 = G1 - dG; H2 = H1 - dH; ks2 = ks1 - dks; Ks2 = Ks1 - dKs;

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

dS1 = MeanDynamicsFunction(G1,H1,J2,L1,RE,b,cs,g1,h1,hs,ks1,mu,nu,sSun);
dS2 = MeanDynamicsFunction(G2,H2,J2,L2,RE,b,cs,g2,h2,hs,ks2,mu,nu,sSun);

dim = constants.STTDim;
STTdot1 = getJOrder1(G1,H1,J2,L1,RE,b,cs,g1,h1,hs,ks1,mu,sSun);
STM1 = reshape(S(dim+1:(dim+dim^2)),dim,dim)';
phiDot = STTdot1 * STM1;
phiDotVec1 = reshape(phiDot',dim^2,1);

if constants.STTOrder == 2
    STTdot2 = getJOrder2(G1,H1,J2,L1,RE,b,cs,g1,h1,hs,ks1,mu,sSun);
    phiDotVec1 = [phiDotVec1; Compute2ndOrderSTTdot(S((dim+dim^2+1):584), STTdot2, dim, STTdot1, STM1)];
end

STTdot1 = getJOrder1(G2,H2,J2,L2,RE,b,cs,g2,h2,hs,ks2,mu,sSun);
STM1 = reshape(S(593:656),dim,dim)';
phiDot = STTdot1 * STM1;
phiDotVec2 = reshape(phiDot',dim^2,1);

if constants.STTOrder == 2
    STTdot2 = getJOrder2(G2,H2,J2,L2,RE,b,cs,g2,h2,hs,ks2,mu,sSun);
    phiDotVec2 = [phiDotVec2; Compute2ndOrderSTTdot(S(657:end), STTdot2, dim, STTdot1, STM1)];
end

dSSTTpropagation = [dS1;phiDotVec1;dS2;phiDotVec2];

end
