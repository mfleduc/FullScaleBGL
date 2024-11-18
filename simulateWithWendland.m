%%%%%%%%
cd('/glade/work/mleduc/GP emulation/needlet fbgl');
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
% load('FSBGL results 3 levels time 1 lambda 1.0000 covmodel mix .mat');
% load('FSBGL results 3 levels time 1 lambda 1.0000 covmodel GCW.mat');
load('FSBGL results 4 levels time 1 lambda 0.1000 covmodel 3 MixOfWendlandCov.mat');
% load('FSBGL results 4 levels time 1 lambda 1.0000 covmodel 3 MixOfWendlandCov.mat');
% load('FSBGL Nugget first/FSBGL results 3 levels time 1 lambda 0.2500 covmodel 3 MixOfWendlandCov nugg first.mat')
load(sprintf('needlets degree resolution j=5.mat'), 'A');

load('mean of beta algorithm residual transform alt final.mat' ,'coeffs', ...
    'a00','c00', 'testdata','finalfit','LAT', 'LON', 'lambda' );
[XH, YH] = sph2hammer(pi/180*(LON-180),pi/180*LAT);
% lambdas = [10:-2:2,0.5];
Sigmacoeffs = inv(Qest{2});
csim = mvnrnd([ones(size(Sigmacoeffs,1 ),1)], Sigmacoeffs, 16);
lowrank = A(:,1:size(Sigmacoeffs,1 )-1)*csim(:,2:end)'+csim(:,1)'/2/sqrt(pi);
params = np(1,:);%[paramrange(1, val(1)),paramrange(2, val(2)), paramrange(3, val(3))];
clear A;
% spatialpart = randn([64800,16]);
% fspart = zeros(size(spatialpart));
[spx,spy,spz] = sph2cart(pi/180*LON, pi/180*LAT,1);
spc = [spx,spy,spz] ;
rows=[];
cols = [];
vals= [];
for ii = 1:length(LAT)
    % d = distance([LAT((ii)),LON((ii))],...
    % [LAT(ii:end), LON(ii:end)]);
    d = real(acosd( spc(ii:end,:)*spc(ii,:)' ));
    kern = WendlandRBF(params(2),d)*params(3)+...
        WendlandRBF(params(4),d)*params(5)+WendlandRBF(params(6),d)*params(7);
    % kern = gentemp([np(2,1:2),-0.1],2*sind(d/2)) + WendlandRBF(np(2,4),d)*np(1,5);
    % d=2*sind(d/2);kern =  gentemp([np(1,1:2),-0.1],d) + gentemp([np(3:4),-0.1],d);
    % kern =  WendlandRBF(params(2),d)*params(3);%+(1-params(3))*WendlandRBF(params(4),d);
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
Sigma = (sparse(rows,cols,vals,64800,64800));
clear rows vals cols
Sigma = (Sigma+Sigma')-(1-10^params(1))*spdiags(diag(Sigma),0,64800,64800);
% Sigma = Sigma + (params(3)+params(5)+10^params(1))*speye(64800) ;
%%%
highres = mvnrnd(zeros(64800,1), Sigma, 16)';
XH = reshape(XH,360,180);
YH = reshape(YH,360,180);
figure;
for ii = 1:3
ax=subplot(2,2,ii); 
% imagesc(flipud(reshape(diag(std(testdata,[],2))*(lowrank(:,ii)+highres(:,ii)) ,360,180)'));colorbar;
pcolor(XH,YH,reshape(diag(std(testdata,[],2))*(lowrank(:,ii)+highres(:,ii)),360,180))
caxis([-5,5]);title('Simulation');shading flat;colorbar;ax.Visible='off';
end
ax=subplot(2,2,4);
% imagesc(flipud(reshape(testdata(:,2,1)-mean(testdata,2), 360,180)'));colorbar;caxis([-5,5])
pcolor(XH,YH,reshape(testdata(:,2,1)-mean(testdata,2) , 360,180));shading flat
title('Data');caxis([-5,5]);colorbar;ax.Visible='off';
 sgtitle( sprintf('3 Wendlands, FS-BGL, 4 levels, decorr length = %.2f deg',params(2)) )
% sgtitle( sprintf('Mix of GC, FS-BGL, 4 levels, decorr length = %.2f deg',4*asind(0.5*params(2))) )
savefig('FSBGL results/simulation 3 levels 3 wendlands nugget first.fig')
%figure;pcolor(XH,YH,(reshape(diag(std(testdata,[],2))*(lowrank(:,ii)+highres(:,ii)) ,360,180)'))