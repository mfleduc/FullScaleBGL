%%%%%%%%
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/NeedMat/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/Spherical-Harmonic-Transform/') );
addpath( genpath('/glade/work/mleduc/Ionosphere_inversion/spherepts/') );
addpath( genpath('/glade/work/mleduc/BigQUIC_release/'))
addpath( genpath('/glade/work/mleduc/MultivariateBasisGraphicalLasso'))
addpath( genpath('/glade/work/mleduc/FMGL-0/'))
addpath( genpath('/glade/work/mleduc/GP emulation/needlet fbgl/NeedletBGL'));
addpath( genpath('/glade/work/mleduc/S2-Sampling-Toolbox'));


% load('FSBGL results 3 levels time 1 lambda 1.0000 covmodel WendlandCov.mat')
% load('FSBGL results 4 levels time 1 lambda 1.0000 covmodel gentemp.mat');
% load('FSBGL results 3 levels time 1 lambda 0.1000 covmodel MixOfWendlandCov.mat');
% load(sprintf('needlets degree resolution j=5.mat'), 'A');
% 
load('mean of beta algorithm residual transform alt final.mat', ...
    'a00','c00', 'testdata' );
load('Fit D results Mix ofGC 4 lvls.mat');
residuals = testdata-finalfit;

[XH, YH] = sph2hammer(pi/180*(LON-180),pi/180*LAT);
XH = reshape(XH,360,180);
YH = reshape(YH,360,180);
% lambdas = [10:-2:2,0.5];
Sigmacoeffs = inv(inputs.Qest{1});
csim = mvnrnd([ones(size(Sigmacoeffs,1 ),1)], Sigmacoeffs, 16);
stdevs = [std(a00(:,1));std(inputs.coeffs,[],2)];
csim = diag(stdevs)*csim' ;
lowrank = A(:,1:size(Sigmacoeffs,1 )-1)*csim(2:end,:)+csim(1,:)/2/sqrt(pi);
clear A


params = [paramrange(1, val(1)),paramrange(2, val(2)), paramrange(3, val(3)),paramrange(4,val(4))];

% % spatialpart = randn([64800,16]);
% fspart = zeros(size(spatialpart));
% Sigma = zeros(64800) ; 
[spx,spy,spz] = sph2cart(pi/180*LON, pi/180*LAT,1);
spc = [spx,spy,spz] ;
rows=[];
cols = [];
vals= [];
for ii = 1:length(LAT)
    % d = distance([LAT((ii)),LON((ii))],...
    % [LAT(ii:end), LON(ii:end)]);
    d = real(acosd( spc(ii:end,:)*spc(ii,:)' ));
    d = 2*sind(d/2);
    kern =  gentemp([params(1:2),-0.1],d)+gentemp([params(3:4),-0.1],d);
    % kern =  params(3)*WendlandRBF(params(2),d)+(1-params(3))*WendlandRBF(params(4),d);
    % Sigma(ii,ii:end) =  Sigma(ii,ii:end)+kern' ;
    % if ii == 7000
    %     Sigma = sparse(Sigma);
    % end
    ndcs = find(kern~=0);
    if ii==1
        vals(1:length(ndcs)) = real(kern(ndcs));
        rows(1:length(ndcs)) = ii;
        cols(1:length(ndcs)) = (ii-1)+ndcs;
    else
        vals(end+1:end+length(ndcs)) = real(kern(ndcs));
        rows(end+1:end+length(ndcs)) = ii;
        cols(end+1:end+length(ndcs)) = (ii-1)+ndcs;
    end
end
Sigma = real(sparse(rows,cols,vals,64800,64800));
Sigma = (Sigma+Sigma')-spdiags(diag(Sigma),0,64800,64800);
Sigma = Sigma + (params(1))*speye(64800) ;
%%%
highres = mvnrnd(zeros(64800,1), Sigma, 16)';
residuals = testdata-finalfit;
spatialstdev = std(residuals,[],2);
highres = diag(spatialstdev)*highres;

figure;
for ii = 1:3
subplot(2,2,ii)
% imagesc(flipud(reshape(diag(std(testdata,[],2))*(lowrank(:,ii)+highres(:,ii)) ,360,180)'));colorbar;
pcolor(XH,YH,reshape(lowrank(:,ii)+highres(:,ii)  ,360,180));colorbar;caxis([-6,6])
title('Simulation');shading flat
end
subplot(2,2,4)
% imagesc(flipud(reshape(testdata(:,2,1)-mean(testdata,2), 360,180)'));colorbar;caxis([-5,5])
pcolor(XH,YH,reshape((finalfit(:,4)-mean(finalfit,2))+(residuals(:,4)-mean(residuals,2)), 360,180));shading flat
title('Data');colorbar;caxis([-6,6])
sgtitle( sprintf('3 levels, Mix of Wendland, split model',params(2)) )
savefig(['Split Model results/split model simulation 3 levels Mix Of Wendland.fig'])
%figure;pcolor(XH,YH,(reshape(diag(std(testdata,[],2))*(lowrank(:,ii)+highres(:,ii)) ,360,180)'))