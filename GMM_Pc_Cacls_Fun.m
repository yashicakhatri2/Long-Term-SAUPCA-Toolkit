function Pc = GMM_Pc_Cacls_Fun(so, sf, Pinertiali, R_star)
    C = NaN(3,3);
    xyzhat = [1 0 0; 0 1 0; 0 0 1];
    bhat(1,:) = sf(1:3) / norm(sf(1:3));
    bhat(2,:) = sf(4:6) / norm(sf(4:6));
    bhat(3,:) = cross(bhat(1,:),bhat(2,:)) / norm(cross(bhat(1,:),bhat(2,:)));
    for i = 1:3
        for j = 1:3
            C(i,j) = dot(bhat(i,:),xyzhat(j,:));     
        end
    end
    
    Xnewi = [C zeros(3,3); zeros(3,3) C] * so;
    Pnewi = [C zeros(3,3); zeros(3,3) C] * Pinertiali * [C zeros(3,3); zeros(3,3) C]';

    P_star = Pnewi([1 3],[1 3]);
    ginv = P_star\eye(2);
    R = R_star;
    mu = Xnewi([1,3]);
    n = R; m = -n; n2 = @(x) sqrt(R.^2 - x.^2); m2 = @(x) -1 .* sqrt(R.^2 - x.^2);
    f2 = 1 / (2 * pi * sqrt(det(P_star))) ; tol = 1.0E-10;
    a = ginv(1,1); b = ginv(1,2) ; c = ginv(2,1) ; d = ginv(2,2);

    int1= @(x, z)exp(-.5 * (a .* (x-mu(1)).^2 + d .* (z-mu(2)).^2 + (x-mu(1)).*(z-mu(2)).*(b+c))) ;
    out2 = integral2(int1, m, n , m2, n2,'Abstol',tol) ;
    Pc = f2 * out2;

end