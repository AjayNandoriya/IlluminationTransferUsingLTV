function U = solve_l1(img)
    lambda = 0.1;
    nNeighbors = 16;
    I = log(img+0.1);
    
%     u = Graph_anisoTV_L1_v2_consistent_weights(uint8(I*255/log(2)), lambda, nNeighbors);
    U = Graph_anisoTV_L1_v2_consistent_weights(uint8(img*255), lambda, nNeighbors);
    U = im2double(U);
%     V = I-U;
end