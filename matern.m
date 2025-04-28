function R = matern( dists, nu, range )
R = 2^(1-nu)/gamma(nu)*(dists/range).^nu.*besselk(nu, dists/range);
R(isnan(R)) = 1;
end