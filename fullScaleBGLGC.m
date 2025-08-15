function fullScaleBGLGC(input)
%% Full scale basis graphical lasso
% clear variables; close all;
fprintf('Running Full-Scale Basis Graphical Lasso with %d levels of resolution and lambdas %.4f to %.4f\n'...
    ,input.lvls,input.lambda(1),input.lambda(end))
fprintf('*******GASPARI-COHN COVARIANCE********\n');
%% Example optimization parameters
optparams = struct; 
optparams.algorithm='simanneal';
optparams.ranges = [[-3,3];[0.001,20];[.0001,1]];
optparams.maxiters = 2000; %number of simulated annealing steps
optparams.npts = 501; %Number of lattice points
optparams.x0 = [0,2,0.25];%Initial guess
optparams.pchange = 0.9;
%%
fprintf('********************Fitting residual model**********************\n');

[spx,spy,spz] = sph2cart(pi/180*input.LON(:), pi/180*input.LAT(:),1);
spc=[spx,spy,spz];
d = real(acosd(spc*spc')) ; %d = d - diag(diag(d));
%imaginary values on the diagonal sometimes        
d = 2*sind(d/2); %chordal distance required for gaspari cohn correlation fn
covmodel = @(x)gentemp([x(2:3),-0.1],d,x(1));%+gentemp([x(3:4),-0.1],d);
% covmodel = @(x)gentemp([x(1:2),-0.1], 2*sind(d/2)) + WendlandCov(x(3:end), d, 'cov');
% [val,obj,np ] = CalculateInitialGuess(covmodel , input.Psi ,input.data, optparams);
% fprintf('Done\n');
% alpha=np(1);
fnname = func2str(covmodel );
fnname = fnname(5:end);
ndx = strfind(fnname, '('); 
fnname = fnname(1:ndx(1)-1);
disp(toc)
alpha = 0.24;
np = [0.240098, -0.564, 0.760962, 0.210079];
for l1 = 1:length(input.lambda) %% cross-validate for each lambda
    lambda = input.lambda(l1);
    % Qest = cell(1,1);
    % np = [];
    fprintf('lambda = %.6f, calculating Q|Sigma\n',lambda);
    %% We want to fit a Wendland covariance
    
    [Qest ] = CalculateQ(input.data, input.Psi, input.lambda(l1), ...
        2e-2,covmodel ,np, strcmpi(fnname,'MarkovRandomField') );
    % svdir = 'FSBGL results/temperature/GC/';
    % mkdir(svdir);
     save(sprintf([input.svDir, ...
         'FSBGL results %d levels lambda %.4f covmodel %s.mat'],input.lvls,lambda,fnname ),...
    'Qest', 'np', 'optparams', 'covmodel','lambda','alpha','input','-v7.3')
     disp(toc)
end

end

