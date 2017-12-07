function demo_TVL1

%% parameters
filename = 'cameraman.tif';
lambda = 1;
nNeighbors = 4;     % use 4 or 16 for 2D images

%% main
root=pwd;
ver='v2.3';

f_8bit = imread(filename); 

cd(ver);
u1 = Graph_anisoTV_L1_v2(f_8bit,lambda,nNeighbors,2);   % do not change the last input argument
cd(root);

figure;
subplot(221);image(f_8bit);colormap(gray(256));axis equal tight; title('Original');
subplot(223);image(u1);colormap(gray(256));axis equal tight; title('TV/L^1 Output u');
subplot(224);image((double(f_8bit)-double(u1))+128);colormap(gray(256));axis equal tight; title('TV/L^1 Output v')
