function [nq,kq,epsq]=getVO2nk(lambdaQuery,Tquery,branch)
% getVO2nk  Louisos 2026 VO2 温度依赖光学常数
% lambdaQuery：波长(μm)，可为标量或向量
% Tquery：温度(°C)，范围35~85
% branch：'Heating'或'Cooling'
% 输出：n、k及复介电常数epsilon=(n+ik)^2

persistent H C Tgrid
if isempty(H)
    rootDir=fileparts(mfilename('fullpath'));
    file=fullfile(rootDir,'VO2_Louisos2026_optical_constants_reconstructed.xlsx');
    if ~isfile(file),error('找不到VO2数据库：%s',file);end
    Tgrid=[35 40 45 50 52.5 55 57.5 60 62.5 65 67.5 70 72.5 75 77.5 80 82.5 85];
    H=readBranch(file,'Heating',Tgrid);
    C=readBranch(file,'Cooling',Tgrid);
end

if Tquery<Tgrid(1) || Tquery>Tgrid(end)
    error('温度必须位于35~85 °C，当前不允许外推。');
end

switch lower(strtrim(branch))
    case 'heating', D=H;
    case 'cooling', D=C;
    otherwise, error('branch只能为''Heating''或''Cooling''。');
end

if any(lambdaQuery<D.lambda(1) | lambdaQuery>D.lambda(end))
    error('波长必须位于%.3f~%.3f μm。',D.lambda(1),D.lambda(end));
end

% 先沿温度方向插值
nT=interp1(Tgrid,D.n.',Tquery,'linear').';
kT=interp1(Tgrid,D.k.',Tquery,'linear').';

% 再沿波长方向插值
nq=interp1(D.lambda,nT,lambdaQuery,'linear');
kq=interp1(D.lambda,kT,lambdaQuery,'linear');
epsq=(nq+1i*kq).^2;
end

function D=readBranch(file,sheet,Tgrid)
M=readmatrix(file,'Sheet',sheet);
M=M(isfinite(M(:,1)),:);
M=M(:,1:1+2*numel(Tgrid));
if any(~isfinite(M(:)),'all'),error('%s工作表存在NaN/Inf。',sheet);end
D.lambda=M(:,1);
D.n=M(:,2:2:end);
D.k=M(:,3:2:end);
end