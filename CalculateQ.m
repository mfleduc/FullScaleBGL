function [Qest] = CalculateQ(Y,Phi,  lambda, tolerance,covmodel, noiseparams, precisionflag)
%UNTITLED4 Summary of this function goes here
%   Calculate the precision matrix for the latent process
%TODO: Detailed explanation
% stdevdiag = diag(stdevprofile);
del = tolerance+1;
p = size(Phi, 2); %Number of basis fns
penalty = lambda*(ones(p)-eye(p)) ;
maxiters = 500;
Qs = {} ; 
Qs{1} = eye(p) ;
Sigma = covmodel(noiseparams);
% dbstop if nan
fprintf('Calculating Q|Sigma...');
for kk = 1:maxiters
    fprintf('Linearizing\n')
    Sstar = LinearizeProblem(Y,Phi, Sigma,Qs{kk}, precisionflag );
    fprintf('Solving at iteration %d \n',kk);
    % keyboard
    
    Qs{kk+1} = QUIC(Sstar, penalty, 1e-3);
    if any(isnan(Qs{kk+1}))
       error('NAN encountered in QUIC output, ceasing calclations\n');
    end
     del = norm(Qs{kk+1}-Qs{kk},'fro')/norm(Qs{kk},'fro');
    if del < tolerance
        fprintf('Tolerance reached: del = %.5f \n', del);
        break
    else
        fprintf('Tolerance not reached, del = %.5f\n', del);
    end
end
fprintf('Done\n');

Qest = Qs{kk+1} ;
end
