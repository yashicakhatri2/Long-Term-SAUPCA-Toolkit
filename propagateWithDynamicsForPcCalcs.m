function OUT = propagateWithDynamicsForPcCalcs(COE1,COE2, STMoi, constants, propTime, mu_old1,mu_old2,Q_bar1,Q_bar2)
% Constants
dim = constants.STTDim;
options = constants.options;

% Initial Del Offset calcs
Del1 = [COE_to_Delaunay(COE1,constants.mu,1);0;0];
Del2 = [COE_to_Delaunay(COE2,constants.mu,1);0;0];
if propTime(1) == 0
    Del1 = Del1 + getOffset(Del1,0,constants,0);
    Del2 = Del2 + getOffset(Del2,0,constants,0);
end
Delnorm1 = normalize(Del1,constants.rE,'vec','Del',1);
Delnorm2 = normalize(Del2,constants.rE,'vec','Del',1);
Pc = 0;

% Initiate STT
STMo1 = reshape(STMoi{1}',64,1);
STMo2 = reshape(STMoi{3}',64,1);
PNewCol1 = NaN(64,8);
PNewCol2 = NaN(64,8);
if constants.STTOrder == 2
    for n = 1:dim
        Pmat1(:,:) = STMoi{2}(:,:,n);
        Pmat2(:,:) = STMoi{4}(:,:,n);
        PNewCol1(:,n) = reshape(Pmat1',dim^2,1);
        PNewCol2(:,n) = reshape(Pmat2',dim^2,1);
    end
    STM2Vec1 = reshape(PNewCol1,dim^3,1);
    STM2Vec2 = reshape(PNewCol2,dim^3,1);
    STMo1 = [STMo1;STM2Vec1];
    STMo2 = [STMo2;STM2Vec2];
end

% Propagate State and STT
So = [Delnorm1;STMo1;Delnorm2;STMo2;Pc];
[tf,Sf] = ode113(@(t,S) SDSDynamicsForPcCalcs(t,S,constants, mu_old1, mu_old2,Q_bar1,Q_bar2), propTime, So, options);

OUT = Sf(end,end);

end

