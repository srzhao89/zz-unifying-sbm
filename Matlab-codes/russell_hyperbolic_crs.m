X=[1 1; 1 2; 2 1; 1 1; 1 1; 2 3];
Y=[2 2; 2 2; 2 2; 1 2; 2 1; 1 2];

for k = 1:size(Y,1)
    xk=X(k,:)';
    yk=Y(k,:)';
    [t,fval]=dea_russell_hyperbolic_crs(xk,yk,X,Y);
    russell_hyperbolic_crs(k,1)=fval;
end

round(russell_hyperbolic_crs,3)