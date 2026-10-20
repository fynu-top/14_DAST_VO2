
%% Multilayer optical calculation
clc; clear all; close all; addpath('optical_constants','script');
t1 = clock; 
wave = (400:50:10000)';          % nm
M_mat = nk_data(wave);           % complex index: n - i*k

%% Material IDs:
% 1  MgF2     % 11 Au
% 2  SiO2     % 12 Air
% 3  PDMS     % 13 Ag
% 4  TIO     
% 5  Si3N    
% 6  VO2_I    
% 7  VO2_M     
% 8  HfO2      
% 9  PMMA     
% 10 DAST     

clear layerID
layerID = [2 3 11 3 11];
d_nm    = [90 70 10 100 800];
subID   = 11;

n_layers = M_mat(:,layerID);
E = zeros(size(wave));

for i = 1:numel(wave)
    E(i) = 1 - reflectance_TM(1, M_mat(i,subID), n_layers(i,:),d_nm*1e-9, 0, wave(i)*1e-9);
end


tc(clock, t1);

plot(wave/1000,E,'LineWidth',1.2);
xlabel('Wavelength (\mum)');
ylabel('1-R');
grid on;



