% Demonstrate face image relighting
%
% Li, Qing, Wotao Yin, and Zhigang Deng. "Image-based face illumination 
% transferring using logarithmic total variation models."
% The visual computer 26.1 (2010): 41-49.
%


source_name = 'ss1';
ref_name = 'ref1';
img=  imresize(im2double(imread(['nimg/' source_name '.PNG'])), 1);
shape_pt = load(['nimg/' source_name '.mat']);
shape_pt =double(shape_pt.shape);

imgA=  imresize(im2double(imread(['nimg/' ref_name '.PNG'])), 1);
shape = load(['nimg/' ref_name '.mat']);
shape =double(shape.shape);

results =   relight( img, shape_pt, imgA, shape);
%%
subplot(131);imshow(img);title('Source');
subplot(132);imshow(imgA);title('Reference');
subplot(133);imshow(results.out_img,[]);title('Output');

