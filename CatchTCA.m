function [value,isterminal,direction] = CatchTCA(t,S, constants)
    mu = constants.mu;
    Del1mean = normalize(S(1:8),constants.rE,'vec','Del',0);
    Del2mean = normalize(S(9:16),constants.rE,'vec','Del',0);

    Del1 = Del1mean + getOffset(Del1mean,t*3600,constants,1);
    Del2 = Del2mean + getOffset(Del2mean,t*3600,constants,1);

    S1 = COE_to_Cartesian(COE_to_Delaunay(Del1,mu,0),mu,1);
    S2 = COE_to_Cartesian(COE_to_Delaunay(Del2,mu,0),mu,1);
    
    rvec = S1(1:3) - S2(1:3);
    r = norm(rvec);
    rdot = S1(4:6) - S2(4:6);
    rdtest = dot(rvec / r, rdot);
    
    value = rdtest;
    isterminal = 0;
    direction = 0;
end