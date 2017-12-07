QUICK START:

1. make sure a C compiler has been set up for mex

2. run

>> compile_display

3. if there is no compiler erros, run
>> Graph_anisoTV_L2_v2
>> Graph_anisoTV_L1_v2



FILES:

- compile_display.m

Compile source files into binary code that displays information.
 
- compile_no_display.m

Compile source files into binary code that does not display information.

- compile_display_debug.m

Compile source files in a debug mode into binary code that displays information.

- Graph_anisoTV_L1_v2_consistent_weights.m / Graph_anisoTV_L2_v2_consistent_weights.m

These two files use the new neighbhor weights described in Rice CAAM TR07-09 that are better than those used before in the sense they better approximate the true anisotropic TV.


USAGE: given in the four m-files


CITATION:

D. Goldfarb and W. Yin. Parametric Maximum Flow Algorithmsfor Fast Total Variation Minimization. SIAM Journal on Scientific Computing, 31(5), 3712-3743. [pdf]