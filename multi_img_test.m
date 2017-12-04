
datapath = '\\QCRI-PRECISION-';
NASpath = '\\DS2015XS\Kilimanjaro';
[in_name,b_name,out_name,l] = textread(fullfile(datapath,'/Documents/ajay/CNN/Face/scripts/comb_1508_testing1_different_target.txt'),'%s %s %s %d');
for k=1:length(in_name)
    in_name{k} = fullfile(NASpath,in_name{k}(11:end));
    b_name{k} = fullfile(NASpath,b_name{k}(11:end));
    out_name{k} = fullfile(NASpath,out_name{k}(11:end));
end

[pid,pidname] = textread(fullfile(datapath,'/Documents/ajay/CNN/Face/dataset/pid_list.txt'),'%d %s');

outdir = 'results/in_full';
for k=289:324%1:length(in_name)
% img =  imresize(im2double(imread('/home/qcri/Documents/ajay/CNN/Face/SOA/EdgePreserving2013/code/000.png')),0.5);

if(~exist(b_name{k},'file') || ~exist(out_name{k},'file'))
    disp([b_name{k} ' does not exist']);
    continue;
end

for id=1:length(pid)
    val=strfind(b_name{k},fullfile(pidname{id}));
    if(~isempty(val))
        break;
    end
end

img = im2double(imread(in_name{k}));
img = permute(img,[2 1 3]);
img = imresize(img,0.5);
img(img<0)=0;
img(img>1)=1;

img_pt = load([datapath '/Documents/ajay/CNN/Face/SOA/EdgePreserving2013/code/' sprintf('%03d',pid(id)) '.mat']);
shape_pt = double(img_pt.shape)/2;



imgB = im2double(imread([NASpath '/Dropbox_MIT/MERL_facial/Old_Proc_045-201/s045-041022-02/refl/refl1_' sprintf('%03d',l(k)) '_15.png']));
imgB(imgB<0)=0;
imgB(imgB>1)=1;
imgB = rgb2hsv(imgB);
imgB = permute(imgB(:,:,3),[2 1 3]);
imgB = imresize(imgB,0.5);


ref_pt =load(fullfile(datapath,'/Documents/ajay/CNN/Face/SOA/EdgePreserving2013/code/000.mat'));
shape = double(ref_pt.shape)/2;


%%
tic
results = relight( ...
        img, ...
        shape_pt, ...
        imgB, ...
        shape );


toc

%% compare
imgGT = im2double(imread(out_name{k}));
imgGT = permute(imgGT,[2 1 3]);
imgGT = imresize(imgGT,0.5);
% imgGT = imcrop(imgGT,results.rect);
% img = imcrop(img,results.rect);
% mask = imcrop(results.mask,results.rect);
mask = results.mask;
img_out=results.out_img;
clear img_out2
for c=1:size(imgGT,3)
    imgGT(:,:,c)=imgGT(:,:,c).*mask;
    
    imgt = imgGT(:,:,c);
    h = imhist(imgt(mask));
    imgt = img_out(:,:,c)./max(max(img_out(:,:,c)));
    imgt(mask) = histeq(imgt(mask),h);
    img_out2(:,:,c)=imgt.*mask;
    
end

%% save results

imwrite(imgGT,sprintf('%s/%04d_%03d_%03d_gt.png',outdir,k,id,l(k)));
imwrite(img_out2,sprintf('%s/%04d_%03d_%03d_output.png',outdir,k,id,l(k)));
imwrite(img_out,sprintf('%s/%04d_%03d_%03d_output_wohe.png',outdir,k,id,l(k)));
imwrite(mask,sprintf('%s/%04d_%03d_%03d_mask.png',outdir,k,id,l(k)));
close all;
% subplot(221);imshow(img,[]);title('input');
% subplot(222);imshow(results.out_img,[]);title('output');
% subplot(223);imshow(img_out,[]);title('output histeq');
% subplot(224);imshow(imgGT,[]);title('GT');
% figure(2);imshowpair(imcrop(results.mask,results.rect),img_out)

end

