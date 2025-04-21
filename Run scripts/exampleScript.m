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
%% load in the data
data = load('data/testdata.mat');
input.data = log( data.data(:,1:25)./median((data.data(:,1:25)),2));
input.LAT = data.lat;input.LON = data.lon;
% 25 steps = 5 mins of data. The first 25 steps are what was mostly used in
% the paper, but the whole dataset was used to do comparison with the EOF
% basis
%% Now we can generate the needlets
jMin = 0; %Lowest spatial resolution.
jMax = 3; %Highest spatial resolution. 
B = 2;%Frequency resolution parameter, B=2 is standard
% Need to convert lat/lon to radians and inclination angle for use in
% NeedMat
Psi = get_A(B, jMin, jMax, pi/2 - pi/180*data.lat, pi/180*data.lon, 1000);
Psi = [ones(size(Psi,1),1)/2/sqrt(pi),Psi]; %Adding the constant function
%% Parameters for QUIC/the DC algorithm
lambda = 0.1;%penalty parameter for quic
tol = 0.005;%Tolerance for convergence of quic
svDir = 'results/MOW/';
mkdir(svDir);

%% Generate input struct
input.lambda=lambda;
input.tol=tol;
input.Psi=Psi;
input.lvls = (jMax-jMin)+1;
input.svDir = svDir;  
fullScaleBGLMOW(input);