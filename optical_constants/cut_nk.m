
clear;clc;close all;
% 1. 读取原始数据
filename = '1MgF2.xlsx';
data = readmatrix(filename); 

% 2. 设置抽取间隔（设为10，可将 5666 行缩减至 567 行）
step = 1; 
reduced_data = data(1:step:end, :);

% 3. 保存为新的 Excel 文件
% writematrix(reduced_data, '7VO2_M_reduced.xlsx');

disp(['原始数据行数: ', num2str(size(data, 1))]);
% disp(['减少后数据行数: ', num2str(size(reduced_data, 1))]);