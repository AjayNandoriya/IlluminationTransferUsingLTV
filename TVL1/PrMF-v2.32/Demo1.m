function Demo1
% need MATLAB Image Processing Toolbox
if ~exist('imnoise','file')
    error('need MATLAB Image Processing Toolbox to load image and add noise');
end

disp = true;
filename = 'cameraman.tif';    % 256 x 256
% filename = 'barbara_org.bmp'; % 512 x 512

f_8bit = imnoise(imread(filename),'gaussian',0.01); 
% f_16bit = uint16(double(f_8bit)*257);
lambda1 = 0.05; lambda2=0.05; lambda_length=5;
lambda = (lambda1:lambda2:(lambda1+(lambda_length-1)*lambda2));

starttime = cputime;

% Execution
nNeighbors = 16;
for ii=1:length(lambda)
%     u1 = Graph_anisoTV_L2_v2_consistent_weights(f_8bit, lambda(ii), nNeighbors);
    u1 = Graph_anisoTV_L1_v2_consistent_weights(f_8bit, lambda(ii), nNeighbors);
    if (disp); figure; imshow(u1); title(sprintf('ROF with lambda = %f',lambda(ii))); end
end
time1 = cputime - starttime;

fprintf('Total time of %d calls = %f sec\n', lambda_length, time1);

end