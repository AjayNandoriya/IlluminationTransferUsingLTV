function demo_TVL1_3D

%% parameters
filename = 'cameraman.tif';
lambda = 0.02;
nNeighbors = 6;     % use 6 for 3D images

%% main
root=pwd;
ver='v2.3';

f_8bit = imread(filename);
f_8bit_3D = reshape(repmat(f_8bit(:),32,1),128,128,128);
disp('The testing 3D image has 2097152 voxels');

cd(ver);
u1 = Graph_anisoTV_L1_v2(f_8bit_3D,lambda,nNeighbors,2);   % do not change the last input argument
cd(root);

% figure;
% subplot(221);image(f_8bit);colormap(gray(256));axis equal tight; title('Original');
% subplot(223);image(u1);colormap(gray(256));axis equal tight; title('TV/L^1 Output u');
% subplot(224);image((double(f_8bit)-double(u1))+128);colormap(gray(256));axis equal tight; title('TV/L^1 Output v')
