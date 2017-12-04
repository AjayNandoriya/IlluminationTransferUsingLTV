% Read landmark point from .fp file
%
% Jiansheng Chen, Guangda Su, Jinping He, Shenglan Ben, "Face Image 
% Relighting using Locally Constrained Global Optimization", ECCV'10
%
%01-10-10, Jiansheng Chen

function pt = readfp(fpname)
    
    assert(exist(fpname, 'file')>0);
    
    fid = fopen(fpname);
    pt = fread(fid, [2 105], 'int');
    pt = pt';
    
    % drop some of the points
    droplist = [101, 102, 104, 105, 95, 96, 98, 99, 80, 82, 90, 91, 92, 60, 62, ...
        71, 73, 57, 59, 53, 55, 68, 70, 64, 66, 35, 36, 38, 39, 29, 30, 32, 33, ...
        50, 51, 47, 48, 41, 42, 44, 45, 61, 72];
    droplist = sort(droplist);
    for i = 1:length(droplist)
        pt(droplist(i)-i+1, :) = [];
    end
    
    fclose(fid);

end