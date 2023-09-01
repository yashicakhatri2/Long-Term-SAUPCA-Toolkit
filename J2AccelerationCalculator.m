function aJ2 = J2AccelerationCalculator(State)
rITRF = J2000toITRF(State(1:3),1);
x = rITRF(1);
y = rITRF(2);
z = rITRF(3);
r_ = sqrt(x^2+y^2+z^2);
nmax = 30;
mmax = 30;
s = x/r_; t = y/r_; u = z/r_; sphi = u; cphi = cos(asin(sphi));
r0 = 1; r(1) = s; i0 = 0; i(1) = t;
rho = a/r_; rho0 = mu/r_; rhon(1) = rho * rho0;
for n = 2:nmax
    rhon(n) = rho * rhon(n-1);
end
for m = 2:mmax
    r(m) = s * r(m-1) - t * i(m-1);
    i(m) = s * i(m-1) + t * r(m-1);
end

A(1,1) = sqrt(3) * cphi;
a1 = 0; a2 = 0; a3 = 0; a4 = 0;
for n = 0:nmax
    a11 = 0; a21 = 0; a31 = 0; a41 = 0;
    sum1 = rhon(n+1)/rE;
    for m = 0:mmax
        if n == m
            A(n,m) = cphi * sqrt((2*n+1)/2/n) * A(n-1,n-1);
        else
            A(n,m) = u * sqrt((2*n+1)*(2*n-1)/(n-m)/(n+m)) * A(n-1,m) - sqrt((2*n+1)*(n-m-1)*(n+m-1)/(2*n-3)/(n+m)/(n-m)) * A(n-2,m);
        end
        C = getGeopConstants('C',n,m); % C/S, n, m
        S = getGeopConstants('S',n,m); % C/S, n, m
        cnm1 = sqrt((n-m)*(n+m+1));
        cn1m1 = sqrt((n+m+2)*(n+m+1)/(2*n+3)/(2*n+2));
        D = C * r(m) + S * i(m);
        if m > 1
            E = C * r(m-1) + S * i(m-1);
            F = C * r(m-1) - S * i(m-1);
        elseif m == 0
            E = 0;
            F = 0;
        elseif m == 1
            E = C * r0 + S * i0;
            F = C * r0 - S * i0;
            
        end
        
        a41 = a41 + cn1m1 * A(n+1,m+1) * D;
        a11 = a11 + A(n,m) * m * E;
        a21 = a21 + A(n,m) * m * F;
        a31 = a31 + cnm1 * A(n,m+1) * D;
    end
    a4 = a4 - sum1 * a41;
    a1 = a1 + sum1 * a11;
    a2 = a2 + sum1 * a21;
    a3 = a3 + sum1 * a31;
end



aJ2_ = a4 * rhat + [a1;a2;a3];
aJ2 = J2000toITRF(aJ2_,0);
end