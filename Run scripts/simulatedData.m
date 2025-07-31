%% Run script for the full-scale basis graphical lasso
%Need to add this folder to the path before running!
rng(61619616);
clear variables;close all;clc;
%% Adding paths to required external code
needmatPath = '/home/male7736/Desktop/Research/NeedMat/'; %Path to NeedMat 
sphereptsPath = '/home/male7736/Desktop/Research/spherepts/';%Path to spherepts
shtPath = '/home/male7736/Desktop/Research/Spherical-Harmonic-Transform/';%Path to Spherical-Harmonic-Transform;
quicPath = '/home/male7736/Desktop/Research/MultivariateBasisGraphicalLasso/src/QUIC/';%Path to QUIC;
addpath(genpath(needmatPath));
addpath(genpath(sphereptsPath));
addpath(genpath(shtPath));
addpath(genpath(quicPath));
%% Generate the data
%high res
lat = -85:5:85;
lon = 0:5:359;
[input.LAT,input.LON] = ndgrid(lat,lon);
[spx,spy,spz] = sph2cart(pi/180*input.LON(:), pi/180*input.LAT(:),1);
spc=[spx,spy,spz];
d = real(acosd(spc*spc')) ;
nu = 0.5;
rangeparam = 60;
tausq = 1/16;
sigmasq = 1;
R = sigmasq*matern( d, nu, rangeparam )  ;
R = R.*WendlandRBF( 60, d );
R = R + tausq*eye(numel(input.LAT));
highres = mvnrnd(zeros(1,numel(input.LAT)),R, 100);
%% Now we can generate the needlets
jMin = 0; %Lowest spatial resolution.
jMax = 2; %Highest spatial resolution. 
B = 2;%Frequency resolution parameter, B=2 is standard
% Need to convert lat/lon to radians and inclination angle for use in
% NeedMat
Psi = get_A(B, jMin, jMax, pi/2 - pi/180*input.LAT(:), pi/180*input.LON(:), 1000);
Psi = [ones(size(Psi,1),1)/2/sqrt(pi),Psi]; %Adding the constant function
% Now: generate coefficients
Q = eye( 23 );
Q(2:end,1)=-1/25;
Q(1,2:end)=-1/25;
Q = blkdiagn(Q, 11)/3;


coeffs = mvnrnd(zeros(1,size(Psi,2)),inv(Q), 100);
lowrank = Psi*coeffs';
input.data = highres'+lowrank;
%% Parameters for QUIC/the DC algorithm
lambda = 10.^(-1:1);%penalty parameter for quic
tol = 0.01;%Tolerance for convergence of quic
svDir = 'results/TM/hub_and_spoke';
mkdir(svDir);

%% Generate input struct
input.lambda=lambda;
input.tol=tol;
input.Psi=Psi;
input.lvls = (jMax-jMin)+1;
input.svDir = svDir;  
input.Qtrue=Q;
fullScaleBGLTaperedMatern(input);