% Perform face relighting based on reference images
%
% Jiansheng Chen, Guangda Su, Jinping He, Shenglan Ben, "Face Image 
% Relighting using Locally Constrained Global Optimization", ECCV'10
%
%01-10-10, Jiansheng Chen

% result = img*imgRefA/imgRefB
function results = relight(img, pt, imgRefA, ptRefA)
    results= struct();
    % morphing the reference images
    scaleA = max(size(img,1)/size(imgRefA,1),size(img,2)/size(imgRefA,2));
    imgRefA = imresize(imgRefA, scaleA);    ptRefA = ptRefA*scaleA;
%     scaleB =  max(size(img,1)/size(imgRefB,1),size(img,2)/size(imgRefB,2));
%     imgRefB = imresize(imgRefB,scaleB);     ptRefB = ptRefB*scaleB;
    try
        imgRefA = rgb2gray(imgRefA);
%         imgRefB = rgb2gray(imgRefB);
    catch
    end
    imgRefA = morph(imgRefA, ptRefA, pt); imgRefA = imgRefA(1:size(img,1), 1:size(img,2));
%     imgRefB = morph(imgRefB, ptRefB, pt); imgRefB = imgRefB(1:size(img,1), 1:size(img,2));
    
    % create the mask of convex hull of pt
    [X,Y] = meshgrid(1:size(img, 2), 1:size(img, 1));
    ch = convhulln(pt);
    mask = inpolygon(X,Y,pt(ch(:,1),1),pt(ch(:,1),2));
    margin = 5;
    xpro = sum(mask, 1); ypro = sum(mask, 2);
    left = max(min(find(xpro>0))-margin, 0); right = min(max(find(xpro>0))+margin, size(mask, 2));
    up = max(min(find(ypro>0))-margin, 0); down = min(max(find(ypro>0))+margin, size(mask, 1));
    rect = [left, up, right-left+1, down-up+1];
    
    results.mask= mask;
    results.rect = rect;
    mask = imcrop(mask, rect);
    imgRefA = imcrop(imgRefA, rect).*mask; 
%     imgRefB = imcrop(imgRefB, rect).*mask;
    hsv_image = rgb2hsv(imcrop(img, rect));
%     hsv_image(:,:,1) = hsv_image(:,:,1).*mask;
%     hsv_image(:,:,2) = hsv_image(:,:,2).*mask;
%     hsv_image(:,:,3) = hsv_image(:,:,3).*mask;
    img_crop = hsv2rgb(hsv_image);
    
    % precalculate some of the temp varibles
    win = 3; % size of the shifting window, make it an ODD number!
    color_channel=3;
    [H, W] = size(hsv_image(:,:,color_channel));
    
    addpath(genpath('TVL1\PrMF-v2.32'));
    imgRefA_U = solve_l1(imgRefA,mask);
    img_crop_U = solve_l1(hsv_image(:,:,color_channel),mask);
    
    hsv_image(:,:,color_channel) = imgRefA_U.*hsv_image(:,:,color_channel).*mask./(img_crop_U+0.0001);
    I_out = hsv2rgb(hsv_image);
    
%     I_out = I_out./max(I_out(:));
%     for c=1:size(I_out,3)
%         I_out(:,:,c) = imadjust(I_out(:,:,c));
%     end
    results.out_img = zeros(size(img));
    results.out_img(rect(2)-1 +(1:size(I_out,1)),rect(1)-1+(1:size(I_out,2)),:)=I_out;
%     figure('name', 'Face Image Relighting');
%     subplot(2,2,1); imshow(img_crop); title('Orignal Image');
%     subplot(2,2,2);  imshow(I_out); title('Relighted Image');
%     subplot(2,2,3);  imshow(imgRefA); title('RefA Image');
%     subplot(2,2,4);  imshow(imgRefB); title('RefB Image');
    
end