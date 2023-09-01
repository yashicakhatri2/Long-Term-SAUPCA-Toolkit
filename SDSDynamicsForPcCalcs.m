% SDS Propagation Method
function dSSTTpropagation = SDSDynamicsForPcCalcs(t,S,constants,mu_old1,mu_old2,Q_bar1,Q_bar2)
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
    phiDotVec2 = [phiDotVec2; Compute2ndOrderSTTdot(S(657:1168), STTdot2, dim, STTdot1, STM1)];
end

% Get pc(t)
dim = constants.STTDim;
mu = constants.mu;
rE = constants.rE;
dyn = constants.dynamicsChoice;
if constants.STTOrder == 1
    STMf1{1} = reshape(S(9:72),dim,dim)';
    STMf2{1} = reshape(S(81:144),dim,dim)';
elseif constants.STTOrder == 2
    STMf1{1} = reshape(S(9:72),dim,dim)';
    STMf2{1} = reshape(S(593:656),dim,dim)';
    test1 = reshape(S(73:584),dim^2,dim);
    test2 = reshape(S(657:1168),dim^2,dim);
    for i = 1:dim
        STMf1{2}(:,:,i) = reshape(test1(:,i),dim,dim)';
        STMf2{2}(:,:,i) = reshape(test2(:,i),dim,dim)';
    end
end

Delfnorm1 = S(1:6);
Delfnorm2 = S(585:592);
Delf1 = normalize(Delfnorm1,constants.rE,'vec','Del',0);
Delf2 = normalize(Delfnorm2,constants.rE,'vec','Del',0);
numericSTT1 = STMf1;
numericSTT2 = STMf2;
Del_updatedi1 =  [wrapToPi(Delf1(1:3));Delf1(4:6)];
Del_updatedi2 =  [wrapToPi(Delf2(1:3));Delf2(4:6)];

% Propagating object 1
mu_updated1 = STTcalcs2BP(mu_old1, constants.STTOrder, constants, numericSTT1);
mu_new1(:,1) = normalize(mu_updated1,rE,'vec','Del',0);
New_plot1 = Del_updatedi1 + mu_new1;
if dyn == "SDSplusSP"
    NewOff1 = getOffset([New_plot1;0;0],t*3600,constants,1);
    New_plot1 = New_plot1 + NewOff1(1:6);
end
COE1out = COE_to_Delaunay(New_plot1, mu, 0);
S1 = COE_to_Cartesian_2(COE1out,mu,1);

% Propagating object 2 
mu_updated2 = STTcalcs2BP(mu_old2, constants.STTOrder, constants, numericSTT2);
mu_new2(:,1) = normalize(mu_updated2,rE,'vec','Del',0);
New_plot2 = Del_updatedi2 + mu_new2;
if dyn == "SDSplusSP"
    NewOff2 = getOffset([New_plot2;0;0],t*3600,constants,1);
    New_plot2 = New_plot2 + NewOff2(1:6);
end
COE2out = COE_to_Delaunay(New_plot2, mu, 0);
S2 = COE_to_Cartesian_2(COE2out,mu,1);

sf = S2 - S1; % direction?

[Q_new_1, Q_new_2] = Covariance_Propagation(Q_bar1, Q_bar2, COE1out, COE2out, mu_old1, mu_old2, constants, numericSTT1, numericSTT2); % Propagates both covariances
pct = GMM_Pc_Cacls_Fun(sf, Q_new_1+Q_new_2, constants.R_star,constants.gridPoints); % GMM Probability of Collision Calculator

dSSTTpropagation = [dS1;phiDotVec1;dS2;phiDotVec2;pct];
end
