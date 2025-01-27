% files1 = dir('FSBGL results 3 levels time 1 lambda * covmodel MixOfWendlandCov.mat');
% files1 = fullfile(files1(1).folder,{files1(1:end).name});
% cv1 = 'CV results mix of Wendland 3 levels.mat';
% files2 = dir('FSBGL results 4 levels time 1 lambda * covmodel MixOfWendlandCov.mat');
% files2 = fullfile(files2(1).folder,{files2(1:end).name});
% cv2 = 'CV results mix of Wendland 4 levels.mat';
% CalculateAIC(files1,cv1);
% title('AIC,FSBGL Mix of Wendland')
% CalculateAIC(files2,cv2);
files = dir('/glade/work/mleduc/GP emulation/needlet fbgl/FSBGL Nugget first/MRF/FSBGL results 3 levels*.mat');
% files = dir('/glade/work/mleduc/GP emulation/needlet fbgl/FSBGL Nugget first/nugget cov/BGL results 3 levels*.mat');

% files = dir('/glade/work/mleduc/GP emulation/needlet fbgl/FSBGL Nugget first/PGL/PGL results 3 levels * MixOfWendland.mat');
files = fullfile(files(1).folder,{files(1:end).name});

CalculateAIC(files,[],0)