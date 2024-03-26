function [t, fval] = dea_russell_hyperbolic_crs(xk,yk,X,Y)
[n,p] = size(X); 
q = size(Y,2);
%%%
f = @(t)(sum(t(1:p))+sum(t((p+1):(p+q)).^(-1)))/(p+q);
A = [-diag(xk), zeros(p,q) , X';
     zeros(q,p), diag(yk), -Y'];
b = [zeros(1,p), zeros(1,q)];
lb=[zeros(1, p), ones(1, q), zeros(1, n)];
ub=[ones(1, p),  inf(1, q),  inf(1, n)];
t0=[0.5*ones(1, p),1.5*ones(1, q), ones(1, n)];
[t,fval]=fmincon(f, t0, A, b, [], [], lb, ub);




