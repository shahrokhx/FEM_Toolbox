function [d, rf] = solution(K, R, debc, ebcVals)
dof = length(R);
df = setdiff(1:dof, debc);
Kf = K(df, df);
Rf = R(df) - K(df, debc)*ebcVals;
dfVals = Kf\Rf;
d = zeros(dof,1);
d(debc) = ebcVals;
d(df) = dfVals;
rf = K(debc,:)*d - R(debc);