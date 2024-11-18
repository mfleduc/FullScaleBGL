function [Sstar] = LinearizeProblem(Y,Phi, Sigma,Q, precisionflag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 p = size(Q,1);
if exist('precisionflag','var')==0
    precisionflag = 0;
end
if precisionflag
    SiPhi = Sigma*Phi;
else    
    SiPhi = Sigma\Phi ;
end
projdata = Y'*SiPhi/sqrt(size(Y,2));

PhiSiPhi = Phi'*SiPhi ;

L = chol(Q+PhiSiPhi,'lower');
M = L'\(L\eye(p));M=0.5*(M+M');
% innerpart = Q+PhiSiPhi+projdata'*projdata;
% QpPsPi = inv((Q+PhiSiPhi));
% Sstar = ((Q+PhiSiPhi))\(innerpart/(Q+PhiSiPhi));
Sstar = M+M*(projdata'*projdata)*M;
if any(isnan(Sstar(:)))
    keyboard
end
% Sstar = 0.5*(Sstar+Sstar');
end