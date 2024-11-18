function Corr = gentemp (params,r)    
%
%  Typical usage: c = 3000;  r = 0:1:6000;     
%                 Corr = genspline0(r,c,.5);     (yields Eqn. 4.10 of Gaspari & Cohn (1999)
%                 Corr = genspline0(r,1000,-.1); (yields the dot-dash curve in Fig. 7, GC 1999)
%
%  Code to find the family C0(z,b,c) used to produce Figs. 7&8 of
%  Gaspari & Cohn (1999) "Construction of correlation functions in
%                         two and three dimensions"
%
%  GENSPLINE0 (R,L,B) returns the generalized spline correlation function.
%  This two-parameter function is given by:
%
% C(z) := fi(z/c),   (i-1)c/2 < z <= ic/2,  c = L*sqrt(10/3),  i = 1,2,3,4   (*)% where
%
%   f1(z) = -[48 bz^4 -48 b^2z^4 + 60 b^z^3 - 9 b + 160 b^2z^2 - 40 bz^2
% -40 bz^2 - 66 b^2 -120 b^2z^3 + 20 z^2 + 56 b^2z^5 - 64 bz^5 - 15 z^3 + 24 z^5%   -24 z^4 -3] / [3*norm],   0 <= z <= 1/2,
%
%   f2(z) = [-120 z^2 - 96 z^5 + 32 z^6 - 4 + 48 z + 29 b - 42 b^2 + 80 z^3
%   + 160 b^2z^6 - 192 bz^6 + 60 z^4 + 612 b^2z - 880 bz^3 + 800 b^2z^3
%   - 1080 b^2z^2 + 780 bz^2 - 210 bz - 384 b^2z^5 + 480 bz^5] / [12 z*norm], 1/2 < z <= 1,
%
%   f3(z) = -b [243 - 230 b + 96 bz^6 - 64 z^6 - 720 z^3 + 1620 z^2 - 1134 z
%   + 732 bz - 1320 bz^2 + 240 bz^4
%   + 288 z^5 - 384 bz^5 - 240 z^4 + 800 bz^3]/[12 z*norm],  1 < z <= 3/2,
%
%   f4(z) = 4 b^2 [ -120 z^2 - 16 + 96 z + 40 z^3
%    + 2 z^6 - 12 z^5 + 15 z^4] / [3 z*norm],  3/2 <= z <= 2,    (**)
%
% where   norm = 3 b + 22 b^2 + 1.
%
%       C is the cutoff distance in Eq. (*) .
%       B is the slope parameter in Eq. (**) (dimensionless).
%
%  Greg Gaspari,  March 25, 1996
%  Last Modified:  Nov. 30, 2000    
b=params(3);c=params(2);
norm =  3*b + 22*b*b + 1;
z = r/c;                        % scale to the support [0,2]
Corr = zeros(size(z));          % initialize

i = find(0 <= z & z <= .5);
 Corr(i) = 48*b*z(i).^4 - 48*b*b*z(i).^4 + 60*b*z(i).^3 - 9*b ...
    + 160*b*b*z(i).*z(i) - 40*b*z(i).*z(i) -66*b*b ...
    - 120*b*b*z(i).^3 + 20*z(i).*z(i) + 56*b*b*z(i).^5 ...
    - 64*b*z(i).^5 - 15*z(i).^3 + 24*z(i).^5 - 24*z(i).^4 - 3;
 Corr(i) = - Corr(i)./(3*norm);

i = find(.5 < z & z <= 1.);
 Corr(i) = -120*z(i).*z(i) - 96*z(i).^5 + 32*z(i).^6 - 4 ...
    + 48*z(i) + 29*b - 42*b*b + 80*z(i).^3 + 160*b*b*z(i).^6 ...
    - 192*b*z(i).^6 + 60*z(i).^4 + 612*b*b*z(i) - 880*b*z(i).^3 ...
    + 800*b*b*z(i).^3 - 1080*b*b*z(i).*z(i) + 780*b*z(i).*z(i) ...
    - 210*b*z(i) - 384*b*b*z(i).^5 + 480*b*z(i).^5;
 Corr(i) = Corr(i)./(12*z(i)*norm);

i = find(1. < z & z <= 1.5);
 Corr(i) = - b.* ( 243 - 230*b + 96*b*z(i).^6 - 64*z(i).^6 ...
    - 720*z(i).^3 + 1620*z(i).*z(i) - 1134*z(i) + 732*b*z(i) ...
    - 1320*b*z(i).*z(i) + 240*b*z(i).^4 + 288*z(i).^5 ...
    - 384*b*z(i).^5 - 240*z(i).^4 + 800*b*z(i).^3 );
 Corr(i) =  Corr(i)./(12*z(i)*norm);

i = find(1.5 < z & z <= 2.);
 Corr(i) = 4*b*b.* (-120*z(i).*z(i) - 16 + 96*z(i) ...
     + 40*z(i).^3 + 2*z(i).^6 - 12*z(i).^5 + 15*z(i).^4 );
 Corr(i) = Corr(i)./(3*z(i)*norm);

Corr = params(1)*sparse(Corr);
end

