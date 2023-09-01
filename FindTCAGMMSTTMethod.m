function [Out, GMM_MEAN] = FindTCAGMMSTTMethod(mu_old1, mu_old2, ICState, constants, doptions, Q_bar1, Q_bar2)

% Constants
propTime2 = constants.P_prop - constants.timeWindow;
mu = constants.mu;

% Initial Conditions for TCA iterator
COE1nom = ICState.obj1.COE;
COE2nom = ICState.obj2.COE;
numericSTT1 = constants.initiateSTT;
numericSTT2 = constants.initiateSTT;

%% Propagating to Nom TCA - del
Del1 = propagateWithDynamics(COE1nom, numericSTT1, constants, [0, propTime2/3600], 1);
Del2 = propagateWithDynamics(COE2nom, numericSTT2, constants, [0, propTime2/3600], 1);
numericSTT1 = Del1{2};
numericSTT2 = Del2{2};
propTime1 = propTime2;
propTime2 = constants.P_prop + constants.timeWindow;
COE1 = COE_to_Delaunay(Del1{1}(1:6)',mu,0);
COE2 = COE_to_Delaunay(Del2{1}(1:6)',mu,0);

%% Propagating to Actual TCA
OUT = propagateWithDynamicsForLongTCA(COE1,COE2, {numericSTT1{1},numericSTT1{2},numericSTT2{1},numericSTT2{2}}, constants, [propTime1/3600, propTime2/3600], 1, doptions, 2, mu_old1, mu_old2);
propTime1 = OUT{3} * 3600;
numericSTT1 = OUT{1}{2};
numericSTT2 = OUT{2}{2};
COE1 = COE_to_Delaunay(OUT{1}{1}(1:6)',mu,0);
COE2 = COE_to_Delaunay(OUT{2}{1}(1:6)',mu,0);

Del_updatedi1 =  OUT{1}{1}(1:6)';
Del_updatedi2 =  OUT{2}{1}(1:6)';
t = propTime1/3600;
mu = constants.mu;
rE = constants.rE;

% Correctly offsetting object 1
mu_updated1 = STTcalcs2BP(mu_old1, constants.STTOrder, constants, numericSTT1);
mu_new1(:,1) = normalize(mu_updated1,rE,'vec','Del',0);
New_plot1 = Del_updatedi1 + mu_new1;

NewOff1 = getOffset([New_plot1;0;0],t*3600,constants,1);
New_plot1 = New_plot1 + NewOff1(1:6);

SCOE1 = COE_to_Delaunay(New_plot1,mu,0);
S1 = COE_to_Cartesian(SCOE1,mu,1);

% Correctly offsetting object 2 
mu_updated2 = STTcalcs2BP(mu_old2, constants.STTOrder, constants, numericSTT2);
mu_new2(:,1) = normalize(mu_updated2,rE,'vec','Del',0);
New_plot2 = Del_updatedi2 + mu_new2;

NewOff2 = getOffset([New_plot2;0;0],t*3600,constants,1);
New_plot2 = New_plot2 + NewOff2(1:6);

SCOE2 = COE_to_Delaunay(New_plot2,mu,0);
S2 = COE_to_Cartesian(SCOE2,mu,1);

sf = S2 - S1;
New_plot_COE = SCOE1;
New_plot_COE_2 = SCOE2;
num1 = numericSTT1;
num2 = numericSTT2;

%% Propoagating back to nominal Time
Del1 = propagateWithDynamics(COE1, numericSTT1, constants, [propTime1/3600, constants.P_prop/3600], 1);
Del2 = propagateWithDynamics(COE2, numericSTT2, constants, [propTime1/3600, constants.P_prop/3600], 1);
numericSTT1 = Del1{2};
numericSTT2 = Del2{2};
Del_updatedi1 =  Del1{1}(1:6)';
Del_updatedi2 =  Del2{1}(1:6)';
t = constants.P_prop/3600;

% Correctly offsetting object 1
mu_updated1 = STTcalcs2BP(mu_old1, constants.STTOrder, constants, numericSTT1);
mu_new1(:,1) = normalize(mu_updated1,rE,'vec','Del',0);
New_plot1 = Del_updatedi1 + mu_new1;

NewOff1 = getOffset([New_plot1;0;0],t*3600,constants,1);
New_plot1 = New_plot1 + NewOff1(1:6);

S1 = COE_to_Cartesian(COE_to_Delaunay(New_plot1,mu,0),mu,1);

% Correctly offsetting object 2 
mu_updated2 = STTcalcs2BP(mu_old2, constants.STTOrder, constants, numericSTT2);
mu_new2(:,1) = normalize(mu_updated2,rE,'vec','Del',0);
New_plot2 = Del_updatedi2 + mu_new2;

NewOff2 = getOffset([New_plot2;0;0],t*3600,constants,1);
New_plot2 = New_plot2 + NewOff2(1:6);

S2 = COE_to_Cartesian(COE_to_Delaunay(New_plot2,mu,0),mu,1);

GMM_MEAN = S1;
so = S2 - S1;

%% Covariance Propagation and Pc calcs
[Q_new_1, Q_new_2] = Covariance_Propagation(Q_bar1, Q_bar2, New_plot_COE, New_plot_COE_2, mu_old1, mu_old2, constants, num1, num2); % Propagates both covariances

Out = GMM_Pc_Cacls_Fun(so, sf, Q_new_1+Q_new_2, constants.R_star); % GMM Probability of Collision Calculator               
end
