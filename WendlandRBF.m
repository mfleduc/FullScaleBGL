function phi = WendlandRBF(theta, distmat)
%% Calculating Wendland RBFs on the sphere
d = distmat/theta ;
phi = (1-d).^6.*(35*d.^2+18*d+3)/3 ; 
phi(d>1)=0;
% out=phi;
end