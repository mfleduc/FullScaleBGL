function [output, objvals,pvec] = CalculateInitialGuess(covmodel, Phi, Y, optparams,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
output = [];
% if precisionflag
if ~isempty(varargin)
    precisionflag = varargin{1};
else
    precisionflag=0;
end
objectiveFn = @(x)NoiseObjectiveFirst(x, covmodel,Phi,Y,precisionflag);
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
        if isfield(optparams, 'x0')&&~isempty(optparams.x0)
            for ii = 1:nparams
                [~,x0(ii)] = min(abs(paramranges(ii,:)-optparams.x0(ii)));
            end
            x0=[250,x0];
        else 
            x0=[];
        end
        paramranges = [linspace(0.0001,100,optparams.npts);paramranges];
        [optval, objvals,pvec] = SimulatedAnnealing(objectiveFn, paramranges, optparams.maxiters, optparams.pchange,0.005 ,x0);
        output = [];
        for ii = 1:nparams
            output(ii) = paramranges(ii,optval(ii));
        end
    case 'fmincon'
        [output, objvals] = fmincon(objectiveFn,optparams.x0,[],[],[],[], optparams.ranges(:,1),optparams.ranges(:,2));
end


end