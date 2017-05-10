N = 400;
M = 128;

Comp_DAS_IQ = N^2*M*6;
Comp_SPURS = 2*(16*N+N*(N+log2(sqrt(N))));
Comp_LinInt = N*N*2;


ratio_DAS_2_SPURS = Comp_DAS_IQ/Comp_SPURS
ratio_DAS_2_LinInt = Comp_DAS_IQ/Comp_LinInt