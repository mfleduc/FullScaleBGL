function Sigma = WendlandCov(params, dist, varargin)
%C2 Wendland covariance 
if ~isempty(varargin)
    type = varargin{1};
else
    type = 'cov';
end

theta = params(2);
tausq = 10^params(1);
d = dist/theta; 
phi = (1-d).^6.*(35*d.^2+18*d+3)/3 ;
% phi = (1+6*d).*(1-d).^6;
% phi = 1-1.5*d+0.5*d^3;
phi(d>1)=0;phi = sparse(phi);
% phi = exp(-d);
if strcmpi(type,'corr')
    Sigma = phi*(1-tausq) + tausq*speye(size(phi,1));
else
    Sigma = params(3)*phi ;
    
    if length(varargin) == 2
        stdprof = varargin{2} ;
        stddiag = spdiags(stdprof, 0, size(stdprof,1),size(stdprof,1));
        Sigma = stddiag*Sigma*stddiag ; 
    end
    Sigma = Sigma + tausq*speye(size(phi,1));
    Sigma = 0.5*(Sigma+Sigma');
end