%% Run script for the full-scale basis graphical lasso
%Need to add this folder to the path before running!
function simulatedData(nSamples)
rng(61619616);
disp(cd)
fprintf("%d samples of the simulated field\n", nSamples)
% clear variables;close all;clc;
%% Adding paths to required external code
dirPath = '/home/male7736/Desktop/Research/';
needmatPath = fullfile(dirPath,'NeedMat'); %Path to NeedMat 
sphereptsPath = fullfile(dirPath,'spherepts');%Path to spherepts
shtPath = fullfile(dirPath,'Spherical-Harmonic-Transform');%Path to Spherical-Harmonic-Transform;
quicPath = fullfile(dirPath,'MultivariateBasisGraphicalLasso/src/QUIC/');%Path to QUIC;
addpath(genpath(needmatPath));
addpath(genpath(sphereptsPath));
addpath(genpath(shtPath));
addpath(genpath(quicPath));
%% Generate the data
%high res
% nSamples= 50;
% lat = -85:3:85;
% lon = 0:3:359;
[input.LAT,input.LON] = ndgrid( 0:0.025:1,0:0.025:1  );%ndgrid(lat,lon);
% [spx,spy,spz] = sph2cart(pi/180*input.LON(:), pi/180*input.LAT(:),1);
% spc=[spx,spy,spz];
% d = real(acosd(spc*spc')) ;
d = sqrt(abs(input.LAT(:)-input.LAT(:)').^2+abs(input.LON(:)-input.LON(:)').^2);
nu = 0.5;
rangeparam = 0.3;
tausq = 1/100;
sigmasq = 1;
R = sigmasq*matern( d, nu, 0.025 )  ;
R = R.*WendlandRBF( rangeparam , d );
R = R + tausq*eye(numel(input.LAT));
%% Now we can generate the needlets
jMin = 0; %Lowest spatial resolution.
jMax = 1; %Highest spatial resolution. 
B = 2;%Frequency resolution parameter, B=2 is standard
% % Need to convert lat/lon to radians and inclination angle for use in
% % NeedMat
% Psi = get_A(B, jMin, jMax, pi/2 - pi/180*input.LAT(:), pi/180*input.LON(:), 1000);
% Psi = [ones(size(Psi,1),1)/2/sqrt(pi),Psi]; %Adding the constant function
Psi = [];
for mm = 0:10
    for nn = 0:10
        Psi(:,end+1) = cos( 2*pi*mm*input.LAT(:) ).*cos( 2*pi*nn*input.LON(:) )*0.5;
    end
end

% Q = Q/5;
Q = eye(11)/2 - diag(ones(1,7),4)/2.2 -  diag(ones(1,7),-4)/2.2;
% Q(2:end,1) = -0.22;
% Q(61,:) = 0.1;
% Q(61,61) = 1.25;%Q=Q/2;
% Q= Q+Q';
Q = blkdiagn(Q,11);
% Q = load("results/TM/rand_Q/Q.mat");Q=Q.Q;
for ensnum = 1:50
highres = mvnrnd(zeros(1,numel(input.LAT)),R, nSamples);

coeffs = mvnrnd(zeros(1,size(Psi,2)),inv(Q), nSamples);
lowrank = Psi*coeffs';
input.data = highres'+lowrank;
%% Parameters for QUIC/the DC algorithm
lambda = 10.^(-1:0.5:1);%penalty parameter for quic
tol = 0.005;%Tolerance for convergence of quic
svDir = sprintf('results/TM/smallscale/%d_samples/ens_%d', nSamples, ensnum);
mkdir(svDir);

%% Generate input struct
input.lambda=lambda;
input.tol=tol;
input.Psi=Psi;
input.lvls = (jMax-jMin)+1;
input.svDir = svDir;  
input.Qtrue=Q;
input.nSamples = nSamples;
input.params=[sigmasq,tausq,nu,.15,rangeparam];
fullScaleBGLTaperedMatern(input);
end
end