function [optval,ps,paramvec] = SimulatedAnnealing(fn, paramranges, maxiters, pchange, tol, x0)
%Simulated annealing for minimizing likelihood function
numvars = size(paramranges, 1);
npts = size(paramranges, 2);
if exist('x0', 'var')==0 || isempty(x0)
    x0 = randi(npts,1,numvars );
end
optval = 0;
del = 1+tol;
ps = zeros(maxiters+1,1);
for ii = 1:numvars
    paramvec(ii) = paramranges(ii,x0(ii));
end
ps(1) = fn(paramvec);
maxstuck = 10;
stuckcnt = 0;
for ii = 1:maxiters
  
    delind = (randi(7,[1,numvars])-4);%Choosing a random neighbor
    if all(delind==0)
        delind = (randi(7,[1,numvars])-4);
    end
    x1 = x0+delind;
    if any(x1>npts)
        x1(x1>npts) = npts;
    end
    if any(x1<1)
        x1(x1<1) = 1;
    end
    for kk = 1:numvars
        paramvec1(kk) = paramranges(kk,x1(kk));
    end
    tmp = fn(paramvec1) ;
    pmove = pchange*exp(-sqrt(ii)*100/abs(ps(1))*max(0, (tmp-ps(ii))));
    if rand(1) < pmove
        x0=x1;
        paramvec = paramvec1;
        ps(ii+1) = tmp;
    else
        ps(ii+1) = ps(ii);
    end
    if mod(ii,100)==0
        fprintf(['iteration %d, params [',repmat('%g, ', 1, numel(paramvec)-1),'%g]\n'], ii,paramvec)
    end
    if abs(tmp-ps(ii))/abs(ps(ii)) < tol
        stuckcnt = stuckcnt+1;
        if stuckcnt == maxstuck
            break
        end
    else
        stuckcnt = 0;
    end
end
optval = x0 ;
if ii ~=maxiters
    fprintf(['iteration %d, params [',repmat('%g, ', 1, numel(x0)-1),'%g]\n'], ii,paramvec)
end
ps = ps(1:ii);
end


