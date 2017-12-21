% Perform face relighting based on reference images
%
% Li, Qing, Wotao Yin, and Zhigang Deng. "Image-based face illumination 
% transferring using logarithmic total variation models."
%�The visual computer�26.1 (2010): 41-49.
%


% result = img*L(imgRefA)/L(img)

function results = relight(img, pt, imgRefA, ptRefA)
    results= struct();

    % morphing the reference images
    scaleA = max(size(img,1)/size(imgRefA,1),size(img,2)/size(imgRefA,2));
    imgRefA = imresize(imgRefA, scaleA);    ptRefA = ptRefA*scaleA;
    imgRefA = morph(imgRefA, ptRefA, pt); imgRefA = imgRefA(1:size(img,1), 1:size(img,2),:);
    
    % create the mask of convex hull of pt
    [X,Y] = meshgrid(1:size(img, 2), 1:size(img, 1));
    ch = convhulln(pt);
    mask = inpolygon(X,Y,pt(ch(:,1),1),pt(ch(:,1),2));
    margin = 5;
    xpro = sum(mask, 1); ypro = sum(mask, 2);
    left = max(find(xpro>0, 1 )-margin, 1); right = min(find(xpro>0, 1, 'last' )+margin, size(mask, 2));
    up = max(find(ypro>0, 1 )-margin, 1); down = min(find(ypro>0, 1, 'last' )+margin, size(mask, 1));
    rect = [left, up, right-left+1, down-up+1];
    
    results.mask= mask;
    results.rect = rect;
    mask = imcrop(mask, rect);
    imgRefA = imcrop(imgRefA, rect).*repmat(mask,1,1,size(imgRefA,3)); 
    hsv_image = imcrop(img, rect);
    
    % precalculate some of the temp varibles
    win = 3; % size of the shifting window, make it an ODD number!
    addpath(genpath('TVL1\PrMF-v2.32'));
    for color_channel=1:3
        imgRefA_U(:,:,color_channel) = solve_l1(imgRefA(:,:,color_channel));
        img_crop_U(:,:,color_channel) = solve_l1(hsv_image(:,:,color_channel));

        hsv_image(:,:,color_channel) = imgRefA_U(:,:,color_channel).*(hsv_image(:,:,color_channel)+0.0001).*mask./(img_crop_U(:,:,color_channel)+0.0001);
%         hsv_image(:,:,color_channel) = (imgRefA_U(:,:,color_channel)+ hsv_image(:,:,color_channel)-img_crop_U(:,:,color_channel)).*mask;
    end
    I_out = (hsv_image);
    
    results.out_img = zeros(size(img));
    results.out_img(rect(2)-1 +(1:size(I_out,1)),rect(1)-1+(1:size(I_out,2)),:)=I_out;
    
end