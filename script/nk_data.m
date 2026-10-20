



function M_mat = nk_data(wave)

wave = wave*1e-3;


% addpath('optical_constants');
M = xlsread('1MgF2.xlsx');    M_mat(:,1) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('1MgF2.xlsx');    M_mat(:,2) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('3PDMS.xlsx');    M_mat(:,3) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('4TIO.xlsx');     M_mat(:,4) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('5Si3N4.xlsx');   M_mat(:,5) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('6VO2_I.xlsx');   M_mat(:,6) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('7VO2_M.xlsx');   M_mat(:,7) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('8HfO2.xlsx');    M_mat(:,8) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('9PMMA.xlsx');    M_mat(:,9) = interp1(M(:,1),M(:,2),wave) -  1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('10DAST.xlsx');   M_mat(:,10) = interp1(M(:,1),M(:,2),wave) - 1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('11Au.xlsx');     M_mat(:,11) = interp1(M(:,1),M(:,2),wave) - 1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('12Air.xlsx');    M_mat(:,12) = interp1(M(:,1),M(:,2),wave) - 1i*interp1(M(:,1),M(:,3),wave);
M = xlsread('13Ag.xlsx');    M_mat(:,13) = interp1(M(:,1),M(:,2),wave) - 1i*interp1(M(:,1),M(:,3),wave);


% end







