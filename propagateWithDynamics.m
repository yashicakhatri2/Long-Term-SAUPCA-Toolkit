function OUT = propagateWithDynamics(COE, STMoi, constants, propTime, flagPropagateSTT)
% Constants
dim = constants.STTDim;
options = constants.options;

% Initial Del Offset calcs
Del = [COE_to_Delaunay(COE,constants.mu,1);0;0];
if propTime(1) == 0
    Del = Del + getOffset(Del,0,constants,0);
end
Delnorm = normalize(Del,constants.rE,'vec','Del',1);

if flagPropagateSTT == 1
    % Initiate STT
    STMo = reshape(STMoi{1}',64,1);
    PNewCol = NaN(64,8);
    if constants.STTOrder == 2
        for n = 1:dim
            Pmat(:,:) = STMoi{2}(:,:,n);
            PNewCol(:,n) = reshape(Pmat',dim^2,1);
        end
        STM2Vec = reshape(PNewCol,dim^3,1);
        STMo = [STMo;STM2Vec];
    end
    
    % Propagate State and STT
    So = [Delnorm;STMo];
    [tf, Sf] = ode113(@(t,S) SDSDynamics(t,S,constants), propTime, So, options);
    Sfend = Sf(end,:);

    % Reshape STT for output
    if constants.STTOrder == 1
        STMf{1} = reshape(Sfend(9:end),dim,dim)';
    elseif constants.STTOrder == 2
        STMf{1} = reshape(Sfend(9:72),dim,dim)';
        test1 = reshape(Sfend(73:end),dim^2,dim);
        for i = 1:dim
            STMf{2}(:,:,i) = reshape(test1(:,i),dim,dim)';
        end
    end
    Delfnorm = Sfend(1:6);
    Delf = normalize(Delfnorm,constants.rE,'vec','Del',0);
    OUT{1} = [wrapToPi(Delf(1:3)) Delf(4:6)];
    OUT{2} = STMf;
else
    % Propagate State
    So = [Delnorm];
    [~,Sf] = ode113(@(t,S) SDSDynamicsStateOnly(t,S,constants), propTime, So, options);

    Delfnorm = Sf(end,1:6);
    Delf = normalize(Delfnorm,constants.rE,'vec','Del',0);
    OUT = [wrapToPi(Delf(1:3)) Delf(4:6)];
end

end

