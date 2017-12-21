function compile_display()

root = pwd;
disp('Please wait ... compiling will finish in less than a minute!');

cd('./sources/MultiThread_algorithms');

mex -DMT_TV_L1_1bit -O -output ../../private/mt_TV_L1_1bit mtl1_TVPara_v2.cpp
mex -DMT_TV_L1_1bit -DMATRIX_Lambda -O -output ../../private/mt_TV_L1_1bit_A mtl1_TVPara_v2.cpp
mex -DMT_TV_L1_8bit -O -output ../../private/mt_TV_L1_8bit mtl1_TVPara_v2.cpp
mex -DMT_TV_L1_8bit -DMATRIX_Lambda -O -output ../../private/mt_TV_L1_8bit_A mtl1_TVPara_v2.cpp
mex -DMT_TV_L1_16bit -O -output ../../private/mt_TV_L1_16bit mtl1_TVPara_v2.cpp
mex -DMT_TV_L1_16bit -DMATRIX_Lambda -O -output ../../private/mt_TV_L1_16bit_A mtl1_TVPara_v2.cpp

cd(root);

cd('./sources/DNC_algorithms');

mex -DDNC_TV_L2_1bit -O -output ../../private/dnc_TV_L2_1bit dnc_v2_TVPara.cpp
mex -DDNC_TV_L2_1bit -DMATRIX_Lambda -O -output ../../private/dnc_TV_L2_1bit_A dnc_v2_TVPara.cpp
mex -DDNC_TV_L2_8bit -O -output ../../private/dnc_TV_L2_8bit dnc_v2_TVPara.cpp
mex -DDNC_TV_L2_8bit -DMATRIX_Lambda -O -output ../../private/dnc_TV_L2_8bit_A dnc_v2_TVPara.cpp
mex -DDNC_TV_L2_16bit -O -output ../../private/dnc_TV_L2_16bit dnc_v2_TVPara.cpp
mex -DDNC_TV_L2_16bit -DMATRIX_Lambda -O -output ../../private/dnc_TV_L2_16bit_A dnc_v2_TVPara.cpp

cd(root);

disp('Finish!');