function[val] = NuggetLikelihood(x,Phi_Phi,projdata,trS,n)
 %% Likelihood for the nugget effect model with Q=\alpha I
 % Adapted from code written by Mitchell Krock and used for the fused basis graphical lasso 
l = size(Phi_Phi,1);
m = size(projdata,2);
% 
tausq = 10^x(1);
Phi_Phi_over_nug = Phi_Phi/tausq ;
projdata_over_nug = projdata/(tausq );
cholQ_plus_Phi = chol(x(2)*eye(l)+Phi_Phi_over_nug,'lower');
logdetpart = 2*sum(log(diag(cholQ_plus_Phi))) - l*log(x(2));%ld(Q+P'P/tau^2)-ld(Q)
tracepart  = 1/m*norm( cholQ_plus_Phi\projdata_over_nug,'fro')^2;
constantpart= n*log(tausq ) + trS/(tausq );
val =  logdetpart - tracepart + constantpart;
end

