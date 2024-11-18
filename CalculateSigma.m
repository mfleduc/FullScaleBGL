function [output, objvals] = CalculateSigma(covmodel, Q, Phi, Y, optparams, precisionflag)
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
        if isfield(optparams, 'x0')
            for ii = 1:nparams
                [~,x0(ii)] = min(abs(paramranges(ii,:)-optparams.x0(ii)));
            end
        else 
            x0=[];
        end
        [optval, objvals] = SimulatedAnnealing(objectiveFn, paramranges, optparams.maxiters, optparams.pchange,0.01 ,x0);
        output = [];
        for ii = 1:nparams
            output(ii) = paramranges(ii,optval(ii));
        end
    otherwise 
        error('Only simulated annealing is implemented')
end


end