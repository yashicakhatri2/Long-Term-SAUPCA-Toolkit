function OUT = getICForObj2SDS(COE1, t1, constants, cartOffset)
mu = constants.mu;
options = constants.options;
% constants.dynamicsChoice = "SDSplusSP";
% % Propagate 1 to tf
% Del1 = [COE_to_Delaunay(COE1,mu,1); 0; 0];
% IC1 = COE_to_Delaunay(Del1,mu,0);
% Del = testPropagateWithDynamics(IC1, 0, constants, [0 t1/3600],t1/2/3600, 0, 0);
% SDSS1 = COE_to_Cartesian_2(COE_to_Delaunay(Del{1},mu,0),mu,1);
% MovedSDSS1 = SDSS1 + cartOffset;
% MovedSDSDel1 = COE_to_Delaunay(COE_to_Cartesian_2(MovedSDSS1,mu,0),mu,1);
% IC_2fStar = COE_to_Delaunay(MovedSDSDel1,mu,0)';
% 
% % Propagate 2 back to to
% Del = testPropagateWithDynamics(IC_2fStar, 0, constants, [0 -t1/3600],-t1/2/3600, 0, 0);
% SDSDelS2 = [Del{1}'; 0; 0];
% IC_2oMean = COE_to_Delaunay(Del{1},mu,0);
% IC_2oDel = SDSDelS2 + getOffset(SDSDelS2,-t1/3600,constants,1);
% IC_2o = COE_to_Delaunay(IC_2oDel,mu,0)';
% 
% % Propagate 2 to tf
% Del = testPropagateWithDynamics(IC_2oMean, 0, constants, [0 t1/3600],t1/2/3600, 0, 0); % No offset IC
% SDSS2 = COE_to_Cartesian_2(COE_to_Delaunay(Del{1},mu,0),mu,1);
% 
% Del = propagateWithDynamics(IC_2o,0, constants, [0 t1/3600],t1/2/3600, 0, 0); % Correctly offset IC
% SDSDelS2Cor = COE_to_Cartesian_2(COE_to_Delaunay(Del{1},mu,0),mu,1);
% 
% % Check Difference between final 1 and 2
% norm(SDSS1(1:3) - SDSS2(1:3))
% norm(SDSS1(1:3) - SDSDelS2Cor(1:3))
% OUT = IC_2o;

[~,CartS1] = ode113(@(t,S) dSFull2(t, S, constants), [0 t1], COE_to_Cartesian_2(COE1,mu,1), options);
[~,CartS2o] = ode113(@(t,S) dSFull2(t, S, constants), [t1 0], CartS1(end,:)+cartOffset', options);
[~,CartS2] = ode113(@(t,S) dSFull2(t, S, constants), [0 t1], CartS2o(end,:), options);
norm(CartS1(end,:)-CartS2(end,:))
OUT = COE_to_Cartesian_2(CartS2o(end,:),mu,0);
end

