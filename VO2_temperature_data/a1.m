clear; clc; close all;
%% 用户设置
Tquery=68;                         % 查询温度(°C)，35~85任意温度
branch='Cooling';                  % 'Heating' 或 'Cooling'
lambdaProbe=1.55;                  % 查询波长(μm)
Tplot=[35 40 45 50 60 65 70 82.5 85];
Tfine=35:0.25:85;

%% 数据文件：本m文件与xlsx共同放在VO2_temperature_data文件夹
rootDir=fileparts(mfilename('fullpath'));
if isempty(rootDir), rootDir=pwd; end
xlsxFile=fullfile(rootDir,'VO2_Louisos2026_optical_constants_reconstructed.xlsx');

if ~isfile(xlsxFile)
    error('找不到文件：\n%s\n请把本m文件和XLSX共同放在VO2_temperature_data文件夹。',xlsxFile);
end

Tgrid=[35 40 45 50 52.5 55 57.5 60 62.5 65 67.5 70 72.5 75 77.5 80 82.5 85];

%% 读取 Heating / Cooling
H=loadVO2Sheet(xlsxFile,'Heating',Tgrid);
C=loadVO2Sheet(xlsxFile,'Cooling',Tgrid);

if numel(H.lambda)~=numel(C.lambda) || max(abs(H.lambda-C.lambda))>1e-12
    error('Heating与Cooling波长网格不一致。');
end

%% 校验重建节点
% 40°C = (35°C + 45°C)/2
err40=max([ ...
    max(abs(C.n(:,2)-(C.n(:,1)+C.n(:,3))/2)), ...
    max(abs(C.k(:,2)-(C.k(:,1)+C.k(:,3))/2))]);

% 82.5°C = (80°C + 85°C)/2
err825=max([ ...
    max(abs(C.n(:,17)-(C.n(:,16)+C.n(:,18))/2)), ...
    max(abs(C.k(:,17)-(C.k(:,16)+C.k(:,18))/2))]);

if max(err40,err825)>1e-9
    error(['Cooling重建节点校验失败：40°C误差=%.3e，' ...
           '82.5°C误差=%.3e。'],err40,err825);
end

fprintf('VO2数据读取成功。\n');
fprintf('波长范围：%.3f ~ %.3f μm，共%d个点\n', ...
    H.lambda(1),H.lambda(end),numel(H.lambda));
fprintf(['温度范围：%.1f ~ %.1f °C，共%d个节点\n' ...
    'Cooling重建节点校验通过：40°C、82.5°C\n\n'], ...
    Tgrid(1),Tgrid(end),numel(Tgrid));

%% 选择 Heating / Cooling
switch lower(strtrim(branch))
    case 'heating'
        D=H;
    case 'cooling'
        D=C;
    otherwise
        error('branch只能为''Heating''或''Cooling''。');
end

if Tquery<Tgrid(1) || Tquery>Tgrid(end)
    error('Tquery必须位于%.1f~%.1f °C。',Tgrid(1),Tgrid(end));
end

if lambdaProbe<D.lambda(1) || lambdaProbe>D.lambda(end)
    error('lambdaProbe必须位于%.3f~%.3f μm。', ...
        D.lambda(1),D.lambda(end));
end

%% 任意温度 n、k 和复介电常数
[n_query,k_query]=vo2AtTemperature(D,Tquery);

N_query=n_query+1i*k_query;
eps_query=N_query.^2;

n0=interp1(D.lambda,n_query,lambdaProbe,'linear');
k0=interp1(D.lambda,k_query,lambdaProbe,'linear');
eps0=(n0+1i*k0)^2;

fprintf('当前支路：%s\n',branch);
fprintf('查询温度：%.2f °C\n',Tquery);
fprintf('查询波长：%.4f μm\n',lambdaProbe);
fprintf('n = %.8f\n',n0);
fprintf('k = %.8f\n',k0);
fprintf('epsilon = %.8f %+.8fi\n\n',real(eps0),imag(eps0));

% 后续程序可以直接使用这四个变量
lambda_um=D.lambda;
n_VO2=n_query;
k_VO2=k_query;
epsilon_VO2=eps_query;

%% 指定温度光谱及连续温度响应
[nSpec,kSpec]=vo2AtTemperature(D,Tplot);

[nHfine,kHfine]=vo2AtTemperature(H,Tfine);
[nCfine,kCfine]=vo2AtTemperature(C,Tfine);

nH_probe=interp1(H.lambda,nHfine,lambdaProbe,'linear');
kH_probe=interp1(H.lambda,kHfine,lambdaProbe,'linear');

nC_probe=interp1(C.lambda,nCfine,lambdaProbe,'linear');
kC_probe=interp1(C.lambda,kCfine,lambdaProbe,'linear');

%% 绘图：所有子图放在同一张图
figure('Color','w','Position',[100 80 1200 820]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(D.lambda,nSpec,'LineWidth',1.2);
xlabel('Wavelength (\mum)');
ylabel('n');
title(sprintf('%s: refractive index',branch));
legend(compose('%.1f ^\\circC',Tplot), ...
    'Location','best','FontSize',8);
grid on; box on;
xlim([D.lambda(1) D.lambda(end)]);

nexttile;
plot(D.lambda,kSpec,'LineWidth',1.2);
xlabel('Wavelength (\mum)');
ylabel('k');
title(sprintf('%s: extinction coefficient',branch));
legend(compose('%.1f ^\\circC',Tplot), ...
    'Location','best','FontSize',8);
grid on; box on;
xlim([D.lambda(1) D.lambda(end)]);

nexttile;
plot(Tfine,nH_probe,'-','LineWidth',1.5);
hold on;
plot(Tfine,nC_probe,'--','LineWidth',1.5);
plot(Tgrid,interp1(H.lambda,H.n,lambdaProbe,'linear'), ...
    'o','MarkerSize',4);
plot(Tgrid,interp1(C.lambda,C.n,lambdaProbe,'linear'), ...
    's','MarkerSize',4);
xlabel('Temperature (^{\circ}C)');
ylabel('n');
title(sprintf('n(T) at \\lambda = %.2f \\mum',lambdaProbe));
legend('Heating','Cooling','Heating nodes','Cooling nodes', ...
    'Location','best');
grid on; box on;

nexttile;
plot(Tfine,kH_probe,'-','LineWidth',1.5);
hold on;
plot(Tfine,kC_probe,'--','LineWidth',1.5);
plot(Tgrid,interp1(H.lambda,H.k,lambdaProbe,'linear'), ...
    'o','MarkerSize',4);
plot(Tgrid,interp1(C.lambda,C.k,lambdaProbe,'linear'), ...
    's','MarkerSize',4);
xlabel('Temperature (^{\circ}C)');
ylabel('k');
title(sprintf('k(T) at \\lambda = %.2f \\mum',lambdaProbe));
legend('Heating','Cooling','Heating nodes','Cooling nodes', ...
    'Location','best');
grid on; box on;

title(tl,'Temperature-dependent optical constants of VO_2');

%% 局部函数
function D=loadVO2Sheet(fileName,sheetName,Tgrid)
M=readmatrix(fileName,'Sheet',sheetName);
M=M(isfinite(M(:,1)),:);

ncol=1+2*numel(Tgrid);
if size(M,2)<ncol
    error('%s工作表列数不足：需要%d列，实际%d列。', ...
        sheetName,ncol,size(M,2));
end

M=M(:,1:ncol);

if any(~isfinite(M(:)))
    error('%s工作表存在NaN/Inf。',sheetName);
end

D.lambda=M(:,1);
D.T=Tgrid(:).';
D.n=M(:,2:2:end);
D.k=M(:,3:2:end);
end

function [nq,kq]=vo2AtTemperature(D,Tq)
Tq=Tq(:).';

if any(Tq<D.T(1) | Tq>D.T(end))
    error('温度必须位于%.1f~%.1f °C，不进行外推。', ...
        D.T(1),D.T(end));
end

nq=interp1(D.T,D.n.',Tq,'linear').';
kq=interp1(D.T,D.k.',Tq,'linear').';
end


