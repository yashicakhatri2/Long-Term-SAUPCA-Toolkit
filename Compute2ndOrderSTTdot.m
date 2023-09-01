function phiDotVec = Compute2ndOrderSTTdot(S, STTdot2, dim, AMat1, STM1)
    test1 = reshape(S,dim^2,dim);
    STM2 = NaN(8,8,8);
    for i = 1:dim
        STM2(:,:,i) = reshape(test1(:,i),dim,dim)';
    end
    PNewCol = NaN(64,8);
    phiDot = pagemtimes(AMat1,STM2) + pagemtimes(STTdot2,pagemtimes(STM1,STM1));
    for n = 1:dim
        Pmat(:,:) = phiDot(:,:,n);
        PNewCol(:,n) = reshape(Pmat',dim^2,1);
    end
    phiDotVec = reshape(PNewCol,dim^3,1);
    
end