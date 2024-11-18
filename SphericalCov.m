function Sigma = SphericalCov(params, dist, varargin)
if ~isempty(varargin)
    type = varargin{1};
else
    type = 'cov';
end
theta = params(2);
tausq = params(1);
d = dist/theta; 
% phi = (1-d).^6.*(35*d.^2+18*d+3)/3 ;
phi = 1-1.5*d+0.5*d.^3;
phi(d>1)=0;
% phi = exp(-d);
if strcmpi(type,'corr')
    Sigma = phi*(1-tausq) + tausq*eye(size(phi,1));
else
    Sigma = params(3)*phi + tausq*eye(size(phi,1));
end