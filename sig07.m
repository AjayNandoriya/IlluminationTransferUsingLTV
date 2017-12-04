function I_out = sig07(I_ratio,I_s1,I_r1)
sigma_s = 40;
sigma_d = 0.1;
sigma_i = sigma_d;
winsize = 2;



[H,W,C]=size(I_ratio);
I_out = zeros(size(I_ratio));
for c=1:C
    for i=1:H
        for j=1:W
            ymin = max(i-winsize,1);
            ymax = min(i+winsize,H);
            xmin = max(j-winsize,1);
            xmax = min(j+winsize,W);
            I_ratio_win = I_ratio(ymin:ymax,xmin:xmax,c);
            I_s1_win = I_s1(ymin:ymax,xmin:xmax,c);
            I_r1_win = I_r1(ymin:ymax,xmin:xmax,c);
            ws = fspecial('gaussian',size(I_ratio_win),sigma_s);
            wd = abs(I_s1_win - I_r1_win );
            wd = exp(-wd./sigma_d)+0.00;
            wd = wd./sum(wd(:));

            wi = abs(I_s1_win(winsize+1,winsize+1)-I_r1_win);
            wi = exp(-wi./sigma_i)+0.00;
            wi = wi./sum(wi(:));

            I_out(i,j,c) = sum(sum(wd.*ws.*wi.*I_ratio_win))./(sum(sum(wd.*ws.*wi)));
%             I_out(i,j,c) = sum(sum(I_ratio_win))./(sum(sum(ones(size(I_ratio_win)))));
        end
    end
end

