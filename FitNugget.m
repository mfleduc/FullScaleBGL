function [alpha,tausq] = FitNugget(Phi,data)
%% Fit nugget effect to the data
%%
x0 = [1,1];
lb = [0 0];
ub = [50 50];
n = size(data,1);
m = size(data,2);
tt = Phi'*data;
trS = norm(data,"fro")^2/m;
PhiSPhi = tt*tt'/m;
x = fmincon(@(x) NuggetLikelihood([x],Phi'*Phi,PhiSPhi,trS,n),x0,[],[],[],[],lb,ub);
alpha = x(2);tausq=x(1);
end