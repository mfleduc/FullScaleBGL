function [output,objvals] = CalculateSigma(covmodel, Q, Phi, Y, optparams, precisionflag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
output = [];
% if precisionflag
objectiveFn = @(x)NoiseObjective(x, covmodel,Q,Phi,Y,precisionflag);
% else
%     PQiP = Phi*(Q\Phi');
%     objectiveFn = @(x)NoiseObjectiveQinv( x, covmodel, PQiP, Y ) ;
% end
switch optparams.algorithm
    case'simanneal'
        nparams = size(optparams.ranges, 1);
        paramranges = zeros(nparams, optparams.npts) ;
        
        for ii = 1:nparams
            paramranges(ii,:) = linspace(optparams.ranges(ii,1),optparams.ranges(ii,2), optparams.npts);
        end
        % paramranges(1,:) = 10.^paramranges(1,:);
        if isfield(optparams, 'x0') && ~isempty(optparams.x0)
            for ii = 1:nparams
                [~,x0(ii)] = min(abs(paramranges(ii,:)-optparams.x0(ii)));
            end
        else 
            x0=repmat(ceil(optparams.npts/2),[1,nparams]);
        end
        [optval, objvals] = SimulatedAnnealing(objectiveFn, paramranges, optparams.maxiters, optparams.pchange,0.005 ,x0);
        output = [];
        for ii = 1:nparams
            output(ii) = paramranges(ii,optval(ii));
        end
    case 'fminbnd'
        %Note: Should only use this to fit the Basis Graphical Lasso (or
        %another one-parameter correlation model)
        warning(' Only use the "fminbnd" option if you are exploring a single parameter covariance model like a nugget model')
        [output, objvals] = fminbnd(objectiveFn, optparams.ranges(1),optparams.ranges(2));
    otherwise 
        error('Only simulated annealing ("simanneal") and fminbnd ("fminbnd") are implemented')
end


end