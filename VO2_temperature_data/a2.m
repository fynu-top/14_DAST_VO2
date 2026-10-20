clear; clc; close all;

%% ======================== 设置 ========================
Tgrid=[35 40 45 50 52.5 55 57.5 60 62.5 65 67.5 70 72.5 75 77.5 80 82.5 85];
lambdaCheck=[0.50 1.00 1.55 2.00];     % 热滞回代表波长(μm)
jumpWarnN=1.0;                          % 相邻温度节点最大|Δn|预警阈值
jumpWarnK=1.0;                          % 相邻温度节点最大|Δk|预警阈值
rootDir=fileparts(mfilename('fullpath')); if isempty(rootDir),rootDir=pwd;end
xlsxFile=fullfile(rootDir,'VO2_Louisos2026_optical_constants_reconstructed.xlsx');
if ~isfile(xlsxFile),error('找不到文件：\n%s',xlsxFile);end

%% ======================== 读取数据 ========================
H=loadVO2Sheet(xlsxFile,'Heating',Tgrid);
C=loadVO2Sheet(xlsxFile,'Cooling',Tgrid);
fprintf('\n================ VO2 数据最终验收 ================\n');
fprintf('文件：%s\n',xlsxFile);
fprintf('温度节点：%d 个，%.1f~%.1f °C\n',numel(Tgrid),Tgrid(1),Tgrid(end));
fprintf('Heating：%d 波长点；Cooling：%d 波长点\n',numel(H.lambda),numel(C.lambda));

%% ======================== 基础质量检查 ========================
passGrid=numel(H.lambda)==numel(C.lambda) && max(abs(H.lambda-C.lambda))<1e-12;
passWave=all(diff(H.lambda)>0) && all(diff(C.lambda)>0);
passFinite=all(isfinite([H.n(:);H.k(:);C.n(:);C.k(:)]));
passPositive=all(H.n(:)>0) && all(C.n(:)>0) && all(H.k(:)>=0) && all(C.k(:)>=0);
passSize=size(H.n,2)==18 && size(C.n,2)==18 && size(H.k,2)==18 && size(C.k,2)==18;
fprintf('\n[基础检查]\n');
showCheck('Heating/Cooling 波长网格一致',passGrid);
showCheck('波长严格递增',passWave);
showCheck('18 个温度节点列完整',passSize);
showCheck('数据无 NaN/Inf',passFinite);
showCheck('n>0 且 k>=0',passPositive);

%% ======================== 重建节点检查 ========================
err40n=max(abs(C.n(:,2)-(C.n(:,1)+C.n(:,3))/2));
err40k=max(abs(C.k(:,2)-(C.k(:,1)+C.k(:,3))/2));
err825n=max(abs(C.n(:,17)-(C.n(:,16)+C.n(:,18))/2));
err825k=max(abs(C.k(:,17)-(C.k(:,16)+C.k(:,18))/2));
tol=1e-10;
fprintf('\n[重建节点检查]\n');
fprintf('Cooling 40 °C：   max|Δn|=%.3e，max|Δk|=%.3e  %s\n',err40n,err40k,status(max(err40n,err40k)<tol));
fprintf('Cooling 82.5 °C： max|Δn|=%.3e，max|Δk|=%.3e  %s\n',err825n,err825k,status(max(err825n,err825k)<tol));

%% ======================== 介电常数 ========================
epsH=(H.n+1i*H.k).^2; epsC=(C.n+1i*C.k).^2;
epsHr=real(epsH); epsHi=imag(epsH);
epsCr=real(epsC); epsCi=imag(epsC);

%% ======================== 相邻温度跳变扫描 ========================
dHn=max(abs(diff(H.n,1,2)),[],1); dHk=max(abs(diff(H.k,1,2)),[],1);
dCn=max(abs(diff(C.n,1,2)),[],1); dCk=max(abs(diff(C.k,1,2)),[],1);
fprintf('\n[相邻温度节点最大变化，全波段扫描]\n');
fprintf('    区间(°C)        H:Δn     H:Δk     C:Δn     C:Δk\n');
for j=1:numel(Tgrid)-1
    flag=dHn(j)>jumpWarnN || dHk(j)>jumpWarnK || dCn(j)>jumpWarnN || dCk(j)>jumpWarnK;
    fprintf('%5.1f→%-5.1f       %6.3f    %6.3f    %6.3f    %6.3f%s\n',...
        Tgrid(j),Tgrid(j+1),dHn(j),dHk(j),dCn(j),dCk(j),ternary(flag,'  *',''));
end
fprintf('* 表示至少一个全波段最大变化超过当前预警阈值；相变区出现较大变化并不等同于数据错误。\n');

%% ======================== 极值检查 ========================
fprintf('\n[全数据范围]\n');
fprintf('Heating：n=%.4f~%.4f，k=%.4f~%.4f，eps''=%.4f~%.4f，eps''''=%.4f~%.4f\n',...
    min(H.n(:)),max(H.n(:)),min(H.k(:)),max(H.k(:)),min(epsHr(:)),max(epsHr(:)),min(epsHi(:)),max(epsHi(:)));
fprintf('Cooling：n=%.4f~%.4f，k=%.4f~%.4f，eps''=%.4f~%.4f，eps''''=%.4f~%.4f\n',...
    min(C.n(:)),max(C.n(:)),min(C.k(:)),max(C.k(:)),min(epsCr(:)),max(epsCr(:)),min(epsCi(:)),max(epsCi(:)));

%% ======================== 总体判定 ========================
basicPass=passGrid && passWave && passFinite && passPositive && passSize;
reconPass=max([err40n err40k err825n err825k])<tol;
fprintf('\n================ 验收结论 ================\n');
if basicPass && reconPass
    fprintf('PASS：基础完整性及 40/82.5 °C 重建校验全部通过。\n');
    fprintf('注意：相邻温度跳变表仅用于发现可疑节点，VO2 相变区本身允许出现快速变化。\n');
else
    fprintf('FAIL：存在基础完整性或重建节点校验问题，请先检查 Excel。\n');
end

%% ======================== 单图 12 子图 ========================
figure('Color','w','Position',[40 40 1550 980]);
tl=tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

plotMap(nexttile,H.lambda,Tgrid,H.n.','Heating: n','n');
plotMap(nexttile,H.lambda,Tgrid,H.k.','Heating: k','k');
plotMap(nexttile,H.lambda,Tgrid,epsHr.','Heating: \epsilon''','\epsilon''');
plotMap(nexttile,H.lambda,Tgrid,epsHi.','Heating: \epsilon''''','\epsilon''''');
plotMap(nexttile,C.lambda,Tgrid,C.n.','Cooling: n','n');
plotMap(nexttile,C.lambda,Tgrid,C.k.','Cooling: k','k');
plotMap(nexttile,C.lambda,Tgrid,epsCr.','Cooling: \epsilon''','\epsilon''');
plotMap(nexttile,C.lambda,Tgrid,epsCi.','Cooling: \epsilon''''','\epsilon''''');

for p=1:numel(lambdaCheck)
    ax=nexttile; lam=lambdaCheck(p);
    nH=interp1(H.lambda,H.n,lam,'linear'); nC=interp1(C.lambda,C.n,lam,'linear');
    kH=interp1(H.lambda,H.k,lam,'linear'); kC=interp1(C.lambda,C.k,lam,'linear');
    yyaxis(ax,'left');
    plot(Tgrid,nH,'-o','LineWidth',1.25,'MarkerSize',3); hold on;
    plot(Tgrid,nC,'--s','LineWidth',1.25,'MarkerSize',3); ylabel('n');
    yyaxis(ax,'right');
    plot(Tgrid,kH,'-o','LineWidth',1.25,'MarkerSize',3); hold on;
    plot(Tgrid,kC,'--s','LineWidth',1.25,'MarkerSize',3); ylabel('k');
    xlabel('Temperature (^{\circ}C)'); title(sprintf('\\lambda = %.2f \\mum',lam));
    grid on; box on; xlim([Tgrid(1) Tgrid(end)]);
    if p==1,legend('n Heating','n Cooling','k Heating','k Cooling','Location','best','FontSize',7);end
end
title(tl,'VO_2 temperature-dependent optical constants: final data validation','FontWeight','bold');

%% ======================== 局部函数 ========================
function D=loadVO2Sheet(fileName,sheetName,Tgrid)
M=readmatrix(fileName,'Sheet',sheetName); M=M(isfinite(M(:,1)),:);
ncol=1+2*numel(Tgrid);
if size(M,2)<ncol,error('%s 工作表列数不足：需要 %d 列，实际 %d 列。',sheetName,ncol,size(M,2));end
M=M(:,1:ncol);
if any(~isfinite(M(:)),'all'),error('%s 工作表存在 NaN/Inf。',sheetName);end
D.lambda=M(:,1); D.T=Tgrid(:).'; D.n=M(:,2:2:end); D.k=M(:,3:2:end);
end

function plotMap(ax,lambda,T,Z,ttl,cbl)
imagesc(ax,lambda,T,Z); axis(ax,'xy'); xlabel(ax,'Wavelength (\mum)'); ylabel(ax,'Temperature (^{\circ}C)');
title(ax,ttl); cb=colorbar(ax); ylabel(cb,cbl); box(ax,'on'); xlim(ax,[lambda(1) lambda(end)]);
end

function showCheck(name,tf)
fprintf('%-34s %s\n',name,status(tf));
end

function s=status(tf)
if tf,s='PASS';else,s='FAIL';end
end

function s=ternary(tf,a,b)
if tf,s=a;else,s=b;end
end