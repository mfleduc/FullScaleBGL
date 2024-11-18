function Sigma = HoleEffect(x, d, stdprof, varargin)
%% C = exp(-(h/d1)^2)cos( h/d2 )+tau^2 delta
p = size(d,1);
p2 = size(stdprof,1);

cov1 = cosd(d*x(3));
cov1(d > (270/x(3))) = 0;
Sigma = sparse((exp(-(d/x(2)))).*cov1);
if p2>1
    stdprof = spdiags(stdprof,0,p2,p2 );
    Sigma = x(4)*stdprof*Sigma*stdprof +x(1)*speye(p); 
else
    Sigma = x(4)*Sigma +x(1)*speye(p); 
end

end