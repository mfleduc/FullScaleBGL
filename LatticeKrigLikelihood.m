function val = LatticeKrigLikelihood(params, tr,Phi,PTP,data)
%%
%params = (rho,alhpa1,...,alphan,kappa,log10(tau^2))
tausq = 10^params(end);
n = size(data,1);p=size(data,2);
B = [];
projdata = 1/tausq*data;
kappa = params(length(tr)+2);
for kk = 1:length(tr)
    B = blkdiag(B,params(kk+1)*sparse( GenerateSARFromMesh(tr{kk},kappa) ));
end
Q = 1/params(1)*(B'*B);
Lq = chol(Q);
Lqptp = chol(Q + 1/tausq*PTP,'lower');

trSDi =  tausq*norm(projdata,'fro')^2*(1/p);

trBig = 1/p*norm( Lqptp\(Phi'*projdata),'fro' )^2;

val = -2*sum(log(diag(Lq))) + n*log(tausq);
val = val + 2*sum(log(diag(Lqptp)))+trSDi - trBig;
end