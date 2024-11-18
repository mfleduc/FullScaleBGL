function Sigma = MixOfWendlandCov(params, dist, varargin)
%C2 Wendland covariance 
if ~isempty(varargin)
    type = varargin{1};
else
    type = 'cov';
end
tausq = 10^params(1); %nugget 
theta1 = params(2); %lengthscale, first one
alpha1 = params(3); %Weighting, first one
theta2 = params(4); %length scale, second one
alpha2 = params(5); %length scale, second one
d1 = dist/theta1; 
d2 = (dist/theta2);

phi1 = (1-(dist/theta1)).^6.*(35*(dist/theta1).^2+18*(dist/theta1)+3)/3 ;
phi2 = (1-(dist/theta2)).^6.*(35*(dist/theta2).^2+18*(dist/theta2)+3)/3 ;
% phi = (1+6*d).*(1-d).^6;
% phi = 1-1.5*d+0.5*d^3;
phi1(dist>theta1)=0;phi1 = sparse(phi1);
phi2(dist>theta2)=0;phi2 = sparse(phi2);
try
    theta3=params(6);
    alpha3=params(7);
    phi3 = (1-(dist/theta3)).^6.*(35*(dist/theta3).^2+18*(dist/theta3)+3)/3 ;
    phi3(dist>theta3)=0;phi3 = sparse(phi3);
catch
    alpha3=0;
phi3=0;
end
% phi = exp(-d);
if strcmpi(type,'corr')
    Sigma = phi*(1-tausq) + tausq*speye(size(phi,1));
else
    Sigma = alpha1*phi1+alpha2*phi2 + alpha3*phi3;
    
    if length(varargin) == 2
        stdprof = varargin{2} ;
        stddiag = spdiags(stdprof, 0, size(stdprof,1),size(stdprof,1));
        Sigma = stddiag*Sigma*stddiag ; 
    end
    Sigma = Sigma + tausq*speye(size(phi1,1));
    Sigma = 0.5*(Sigma+Sigma');
end