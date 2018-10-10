
% This file is used to generate golden input data for unit testing.
% DO NOT MODIFY

% these are used to test linspace()
clear

w1 = linspace(-1,1,9);
w2 = linspace(1,-1,9);
w3 = linspace(-1,1,8);

save linspace_neg1_1_9.dat w1;
save linspace_1_neg1_9.dat w2;
save linspace_neg1_1_8.dat w3;

% these are used to test readVectorFromBinFile()
clear

w = linspace(-1,1);
wT = w';

writeToFile("readVectorBin_neg1_1_100.dat", w);
writeToFile("readVectorBin_neg1_1_100T.dat", wT);

% these are used to test readVectorFromTxtFile()
clear

w = linspace(-1,1);
wT = w';

save readVectorTxt_neg1_1_100.dat w;
save readVectorTxt_neg1_1_100T.dat wT;

% these are used to test readVectorFromTxtFile()
clear

for i = 0:4; for j = 0:4
  m(i+1,j+1) = 17/(i+1+j);
endfor; endfor;

save readMatrixTxt_5x5.dat m;

% these are used to test lgwt()
clear

[x10, w10] = lgwt(10,-1,1);
[x100, w100] = lgwt(100,-13,34);

save lgwt_x10.dat x10;
save lgwt_w10.dat w10;
save lgwt_x100.dat x100;
save lgwt_w100.dat w100;

% these are used to test legendre()
clear

[x10, w10] = lgwt(10,-1,1);

v_10_2 = legendre(x10, 2); % second arg is the polynomial degree
v_10_3 = legendre(x10, 3);
v_10_4 = legendre(x10, 4);
v_10_5 = legendre(x10, 5);
v_10_6 = legendre(x10, 6);
v_10_7 = legendre(x10, 7);

save legendre_10_2.dat v_10_2
save legendre_10_3.dat v_10_3
save legendre_10_4.dat v_10_4
save legendre_10_5.dat v_10_5
save legendre_10_6.dat v_10_6
save legendre_10_7.dat v_10_7

[x100, w100] = lgwt(100,-13,34);

v_100_2 = legendre(x100, 2);
v_100_3 = legendre(x100, 3);
v_100_4 = legendre(x100, 4);
v_100_5 = legendre(x100, 5);
v_100_6 = legendre(x100, 6);
v_100_7 = legendre(x100, 7);

save legendre_100_2.dat v_100_2
save legendre_100_3.dat v_100_3
save legendre_100_4.dat v_100_4
save legendre_100_5.dat v_100_5
save legendre_100_6.dat v_100_6
save legendre_100_7.dat v_100_7

% these are used to test dlegendre()
clear

[x10, w10] = lgwt(10,-1,1);

v_10_2 = dlegendre(x10, 2); % second arg is the polynomial degree
v_10_3 = dlegendre(x10, 3);
v_10_4 = dlegendre(x10, 4);
v_10_5 = dlegendre(x10, 5);
v_10_6 = dlegendre(x10, 6);
v_10_7 = dlegendre(x10, 7);

save dlegendre_10_2.dat v_10_2
save dlegendre_10_3.dat v_10_3
save dlegendre_10_4.dat v_10_4
save dlegendre_10_5.dat v_10_5
save dlegendre_10_6.dat v_10_6
save dlegendre_10_7.dat v_10_7

[x100, w100] = lgwt(100,-13,34);

v_100_2 = dlegendre(x100, 2);
v_100_3 = dlegendre(x100, 3);
v_100_4 = dlegendre(x100, 4);
v_100_5 = dlegendre(x100, 5);
v_100_6 = dlegendre(x100, 6);
v_100_7 = dlegendre(x100, 7);

save dlegendre_100_2.dat v_100_2
save dlegendre_100_3.dat v_100_3
save dlegendre_100_4.dat v_100_4
save dlegendre_100_5.dat v_100_5
save dlegendre_100_6.dat v_100_6
save dlegendre_100_7.dat v_100_7

% these are used to test operator_two_scale()
clear

output23 = OperatorTwoScale(2,3);
output34 = OperatorTwoScale(3,4);
output45 = OperatorTwoScale(4,5);
output26 = OperatorTwoScale(2,6);

save operator_two_scale_2_3.dat output23
save operator_two_scale_3_4.dat output34
save operator_two_scale_4_5.dat output45
save operator_two_scale_2_6.dat output26

% FIXME initial_conditions_1D() data is not generated here
