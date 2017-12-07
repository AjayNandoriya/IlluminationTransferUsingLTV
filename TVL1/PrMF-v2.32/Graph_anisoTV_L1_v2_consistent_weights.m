function u=Graph_anisoTV_L1_v2_consistent_weights(F,lambda,nNeighbors,biThread)
% Solve
%    min TV(u) + lambda * ||u-F||_L1
%   OR
%    min TV(u) + sum_i ( lambda_i |u_i - F_i| )
%
% Difference between 
%     Graph_anisoTV_L1_v2
% and
%     Graph_anisoTV_L1_v2_consistent_weights:
% the latter file use the new neighbhor weights described in Rice CAAM TR07-09 that
% are better than those used before in the sense they better approximate the true
% isotropic TV.
%
% F = input 1/8/16-bit image 2D/3D matrix
%                1-bit F: "logical", values: 0, 1
%                8-bit F: "uint8", values: 0, 1, ..., 255
%               16-bit F: "uint16",values: 0, 1, ..., 65535
% lambda = a positive scalar OR a matrix with positive scalars of double type
% nNeighbors =  4:  anisotropic  4 neighbors for 2D images
%            or 8: anisotropic 16 neighbors for 2D images
%            or 16: anisotropic 16 neighbors for 2D images
%            or 5:  isotropic    5 neighbors for 2D binary images, no longer provided
%            or 6:  anisotropic  6 neighbors for 3D images
% biThread: binary number 'xy' where
%              x={0,1} switch for UP thread
%              y={0,1} switch for DN thread
%           So, set
%           xy = 1 for DN thread;
%           xy = 2 for UP thread;
%           xy = 3 for BOTH threads (now turned off for technical reasons)
%
% Data scaling example:
%           To solve
%             min TV(u) + mu ||u-f||_L1, where u and f are matrices with entries in [0,1],
%           using the 16-bit resolution, apply the following steps:
%
%           bit = 16;
%           nNeighbors = 4;
%           scale = 2^bit - 1;
%           F = uint16(scale*f); % map to 16 positive integers
%           lambda = mu;
%           u = Graph_anisoTV_L1_v2_consistent_weights(F,lambda,nNeighbors,biThread);
%           u = double(u)/scale; % map to [0,1]
%
% Author: Wotao Yin, wotao.yin@rice.edu, 2010

%% test
if ~(exist('F','var') || exist('lambda','var')...
   ||exist('nNeighbors','var') || exist('biThread','var'))
    warndlg('No input is specified. Running a tiny demo.');
    F=logical([0 1 1; 1 0 0]);
    lambda = 1.5;
    nNeighbors = 4;
    biThread = 2;
end

if ~exist('biThread','var'); biThread=2; end

%% Graph neighborhood topology
    %% type information
    % each row correspondes to an arc
    % 1st element: 0 = incoming; 1 = outgoing;
    %              2 = related arcs;
    %              3 = incoming s arc;
    %              4 = outgoing t arc;
    %              5 = related incoming s arc;
    %              6 = ralated outgoing t arc;
    %              Now: only 1, 3, 4 are allowed
    % 2nd element: outgoing arc capacity (except for s-arc, that's the incoming capacity)
    % (3rd,4th) elements: (row col) offset to the partner node if not terminal arc
ndim = ndims(F);
if (any(size(F)==1)); error('There is a trivial dimension (size=1) in F. Please correct.'); end

if (nNeighbors==4 && ndim==2)
    w4=pi/4;
    % 2D: 4-point neighbors
    type = [
        1   w4         0   1 ;
        1   w4         1   0 ;
        1   w4         0  -1 ;
        1   w4        -1   0 ]';

elseif (nNeighbors==8 && ndim==2)
    w4=pi/8; w8=sqrt(2)*pi/16;
    % 2D: 8-point neighbors
    type = [
        1   w4         0   1 ;      
        1   w8         1   1 ;
        1   w4         1   0 ;
        1   w8         1  -1 ;
        1   w4         0  -1 ;
        1   w8        -1  -1 ;
        1   w4        -1   0 ;
        1   w8        -1   1 ]';

elseif (nNeighbors==16 && ndim==2)
    w4=atan(0.5)/2; w8=(pi/4-atan(0.5))/(2*sqrt(2)); w16=(atan(0.5)+(pi/4-atan(0.5)))/2/(2*sqrt(5));
    % 2D: 16-point neighbor
    type = [
        1   w4          0   1;
        1   w4          1   0;
        1   w4          0  -1;
        1   w4         -1   0;
        1   w8          1   1;
        1   w8          1  -1;
        1   w8         -1  -1;
        1   w8         -1   1;
        1   w16         1   2;
        1   w16         2   1;
        1   w16         2  -1;
        1   w16         1  -2;
        1   w16        -1  -2;
        1   w16        -2  -1;
        1   w16        -2   1;
        1   w16        -1   2]';
    
elseif (nNeighbors==6 && ndim==3)
    w6 = (pi/2)^2/pi;
    % 3D: 6-point neighbor
    type = [
        1   w6           0   0   1;
        1   w6           0   1   0;
        1   w6           0   0  -1;
        1   w6           0  -1   0;
        1   w6           1   0   0;
        1   w6          -1   0   0]';
else
    error('Unsupported dimension or neighbor type.');
end
        
    
%% call parametric max-flow
if (isa(F, 'logical'))
    if (isscalar(lambda))
        u = mt_TV_L1_1bit(F,lambda,type,biThread);
    else
        u = mt_TV_L1_1bit_A(F,lambda,type,biThread);
    end
elseif (isa(F, 'uint8'))
    if (isscalar(lambda))
        u = mt_TV_L1_8bit(F,lambda,type,biThread);
    else
        u = mt_TV_L1_8bit_A(F,lambda,type,biThread);
    end
elseif (isa(F, 'uint16'))
    if (isscalar(lambda))
        u = mt_TV_L1_16bit(F,lambda,type,biThread);
    else
        u = mt_TV_L1_16bit_A(F,lambda,type,biThread);
    end
else
    error('This program only supports inputs of type uint8 or uint16.');
end

%% Show final energy for verification
%energy_l2TVL1(F,u,lambda);