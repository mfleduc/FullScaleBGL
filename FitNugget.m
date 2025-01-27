function [alpha,logtausq] = FitNugget(Phi,data)
%% Fit nugget effect to the data
%%
x0 = [-4,1.5];
lb = [-10 0];
ub = [1 1500];
n = size(data,1);
m = size(data,2);
tt = Phi'*data;
trS = norm(data,"fro")^2/m;
[x,val] = fmincon(@(x) NuggetLikelihood(x,Phi'*Phi,tt,trS,n),x0,[],[],[],[],lb,ub);
alpha = x(2);logtausq=x(1);
end