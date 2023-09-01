function J = Partials_PQW_COE(COE, mu, t, check)

a = COE(1); e = COE(2); i = COE(3); O = COE(4); w = COE(5); M_t = COE(6) + sqrt(mu/a^3) * t;
P = [1; 0; 0]; 
Q = [0; 1; 0]; 
R = [0; 0; 1];

TA = MAtoTA(M_t,e); 
E = 2 * atan(sqrt((1-e)/(1+e)) * tan(TA/2));
r = a * (1 - e * cos(E));
n = sqrt(mu/a^3);

X = a * (cos(E) - e); 
Y = a * sqrt(1-e^2) * sin(E);
Xdot = -n * a^2 * sin(E) / r; 
Ydot = n * a^2 * sqrt(1 - e^2) * cos(E) / r;
x = P * X + Q * Y; 
xdot = P * Xdot + Q * Ydot;

L = a^2  * (e * cos(E) - 1 - sin(E)^2) / r; 
M = a^2 * sin(E) * (cos(E) - e) / r / sqrt(1-e^2);
Ldot = n * a^4 * (e - 2 * cos(E) + e * cos(E)^2) * sin(E) / r^3;
Mdot = n * a^4 * (e^2 - 1 - e * cos(E) + 2 * cos(E)^2 - e * cos(E)^3) / r^3 / sqrt(1 - e^2);

dxda = 1/a * (x - 3/2*xdot*t);
dxde = L * P + M * Q;
dxdM = xdot / n;
dxdi = (X * sin(w) + Y * cos(w)) * R;
dxdw = Q * X - P * Y;
dxdO = [-x(2); x(1); 0];
dxdotda = -1/2/a * (xdot - 3 * mu * x * t / r^3);
dxdotde = Ldot * P + Mdot * Q;
dxdotdM = -n * (a/r)^3 * x;
dxdotdi = (Xdot * sin(w) + Ydot * cos(w)) * R;
dxdotdw = Q * Xdot - P * Ydot;
dxdotdO = [-xdot(2); xdot(1); 0];

dadx = 2 * a^2 / r^3 * x;
dedx = sqrt(1-e^2)/n/a^2/e*(Q * Xdot - P * Ydot + n * sqrt(1-e^2) * (a/r)^3 * x);
didx = ((P * Ydot - Q * Xdot) * cos(i) + dxdotdO) / (n * a^2 * sqrt(1-e^2) * sin(i));
dMdx = 1/n/a^2 * (-xdot + 3*mu*x*t/r^3 + (1-e^2)/e*(Ldot*P + Mdot*Q));
dwdx = 1/n/a^2 * (-sqrt(1-e^2)/e*(Ldot*P + Mdot*Q) + cot(i)/sqrt(1-e^2)*dxdotdi);
dOdx = -dxdotdi/(n*a^2*sqrt(1-e^2)*sin(i));
dadxdot = 2 * xdot / n^2 / a;
dedxdot = sqrt(1-e^2)/n/a^2/e*(P*Y - Q*X + sqrt(1-e^2)*xdot/n);
didxdot = -((P*Y - Q*X) * cos(i) + dxdO) / (n*a^2*sqrt(1-e^2)*sin(i));
dMdxdot = 1/n/a^2*(-2*x + 3*xdot*t - (1-e^2)/e*(L*P + M*Q));
dwdxdot = -1/n/a^2*(-sqrt(1-e^2)/e*(L*P + M*Q) + cot(i)*dxdi/sqrt(1-e^2));
dOdxdot = dxdi/(n*a^2*sqrt(1-e^2)*sin(i));



if check == 1
    J = [dxda dxde dxdi dxdO dxdw dxdM;...
        dxdotda dxdotde dxdotdi dxdotdO dxdotdw dxdotdM];
else
    J = [dadx' dadxdot';...
        dedx' dedxdot';...
        didx' didxdot';...
        dOdx' dOdxdot';...
        dwdx' dwdxdot';...
        dMdx' dMdxdot'];
        
end

end