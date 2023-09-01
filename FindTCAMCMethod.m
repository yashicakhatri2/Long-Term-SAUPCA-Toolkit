function [OUT1, OUT2] = FindTCAMCMethod(COE1nom, COE2nom, P_nom, constants)
mu = constants.mu;
propTime1 = 0;
COE1 = COE1nom;
COE2 = COE2nom;

if constants.testFig2 == 1
    propTime2 = P_nom;
    Del1mean = [propagateWithDynamics(COE1, 0, constants, [propTime1/3600, propTime2/3600], 0),0,0];
    Del2mean = [propagateWithDynamics(COE2, 0, constants, [propTime1/3600, propTime2/3600],  0),0,0];
    t = propTime2;
    Del1 = Del1mean' + getOffset(Del1mean,t*3600,constants,1);
    Del2 = Del2mean' + getOffset(Del2mean,t*3600,constants,1);

    OUT1 = COE_to_Cartesian(COE_to_Delaunay(Del1,mu,0),mu,1);
    OUT2 = COE_to_Cartesian(COE_to_Delaunay(Del2,mu,0),mu,1);
else
    propTime2 = P_nom - constants.timeWindow;
    EventFun = @(t,S) CatchTCA(t,S,constants);
    doptions = odeset('RelTol', 1E-13,'AbsTol',1E-13,'Events',EventFun);
    
    
    Del1 = propagateWithDynamics(COE1, 0, constants, [propTime1/3600, propTime2/3600], 0);
    Del2 = propagateWithDynamics(COE2, 0, constants, [propTime1/3600, propTime2/3600],  0);
    COE1 = COE_to_Delaunay(Del1,mu,0);
    COE2 = COE_to_Delaunay(Del2,mu,0);
    
    propTime1 = propTime2;
    propTime2 = P_nom + constants.timeWindow;
    
    OUT = propagateWithDynamicsForLongTCA(COE1, COE2, 0, constants, [propTime1/3600, propTime2/3600], 0, doptions, 1, 1, 1);
    
    OUT1 = OUT{1};
    OUT2 = OUT{2};
end

end
