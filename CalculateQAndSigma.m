function [Qest, noisemodel] = CalculateQAndSigma(Y,Phi,  lambda, tolerance,covmodel, optparams, precisionflag,nuggresults)
%UNTITLED4 Summary of this function goes here
%   Use if you are directly modeling Sigma, not Sigma^{-1}
%TODO: Detailed explanation
% stdevdiag = diag(stdevprofile);
del = tolerance+1;
p = size(Phi, 2); %Number of basis fns
penalty = lambda*(ones(p)-eye(p)) ;
maxiters = 100;
Qs = {} ; 
if exist('nuggresults','var')==1
    Sigma = nuggresults*speye(size(Phi,1));
    Qs{1} = speye(p);
else
    Qs{1} = (1+0.5^2)*eye(p) ;
    
    np0 = optparams.x0;
    np = zeros(maxiters+1, length(np0));
    np(1,:) = np0;
    Sigma = covmodel(np0);
end
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
if ~isempty(optparams.algorithm)
    fprintf('Calculating Sigma|Q...\n');
    np = CalculateSigma(covmodel, Qs{kk+1}, Phi, Y, optparams, precisionflag) ;
    fprintf('Done\n');
else
    np = optparams.x0;
end
Qest = Qs{kk+1} ;noisemodel = np ;
end
