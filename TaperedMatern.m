function R = TaperedMatern(d,x)%Wendland tapered matern correlation function
R = speye( size(d) )*10^x(1);
R = R + x(2)*sparse(double(WendlandRBF(x(3),d).*matern(d,x(4),x(5))));
end