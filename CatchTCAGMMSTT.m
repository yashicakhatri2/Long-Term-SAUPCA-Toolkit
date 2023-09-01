function [value,isterminal,direction] = CatchTCAGMMSTT(t,S, constants,mu_old1,mu_old2)
    
    dim = constants.STTDim;
    mu = constants.mu;
    rE = constants.rE;
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
    Out1 = [wrapToPi(Delf1(1:3)) Delf1(4:6)];
    Out2 = [wrapToPi(Delf2(1:3)) Delf2(4:6)];
    numericSTT1 = STMf1;
    numericSTT2 = STMf2;
    Del_updatedi1 =  Out1(1:6)';
    Del_updatedi2 =  Out2(1:6)';
    
    % Propagating object 1
    mu_updated1 = STTcalcs2BP(mu_old1, constants.STTOrder, constants, numericSTT1);
    mu_new1(:,1) = normalize(mu_updated1,rE,'vec','Del',0);
    New_plot1 = Del_updatedi1 + mu_new1;

    NewOff1 = getOffset([New_plot1;0;0],t*3600,constants,1);
    New_plot1 = New_plot1 + NewOff1(1:6);

    S1 = COE_to_Cartesian(COE_to_Delaunay(New_plot1,mu,0),mu,1);
    
    % Propagating object 2 
    mu_updated2 = STTcalcs2BP(mu_old2, constants.STTOrder, constants, numericSTT2);
    mu_new2(:,1) = normalize(mu_updated2,rE,'vec','Del',0);
    New_plot2 = Del_updatedi2 + mu_new2;
    
    NewOff2 = getOffset([New_plot2;0;0],t*3600,constants,1);
    New_plot2 = New_plot2 + NewOff2(1:6);
    
    S2 = COE_to_Cartesian(COE_to_Delaunay(New_plot2,mu,0),mu,1);
    
    rvec = S1(1:3) - S2(1:3);
    r = norm(rvec);
    rdot = S1(4:6) - S2(4:6);
    rdtest = dot(rvec / r, rdot);
    
    value = rdtest;
    isterminal = 0;
    direction = 0;
end
