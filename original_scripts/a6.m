% exportgraphics(gcf,'DAST_robustness.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

FS1=20;FS2=26;FN='Times New Roman';LW=1.5;FigWidth=1900;FigHeight=880;
set(groot,'defaultAxesFontName',FN,'defaultTextFontName',FN,'defaultLegendFontName',FN,...
    'defaultAxesFontSize',FS1,'defaultTextFontSize',FS1,'defaultLegendFontSize',FS1,...
    'defaultAxesFontWeight','normal','defaultTextFontWeight','normal',...
    'defaultAxesLabelFontSizeMultiplier',1,'defaultAxesTitleFontSizeMultiplier',1,...
    'defaultTextInterpreter','tex','defaultLegendInterpreter','tex',...
    'defaultAxesTickLabelInterpreter','tex');

%% Paths and materials
root=fileparts(mfilename('fullpath'));if isempty(root),root=pwd;end
project=fileparts(root);
P={root,project,fullfile(root,'optical_constants'),fullfile(root,'script'),...
    fullfile(root,'Atmospheric data'),fullfile(project,'optical_constants'),...
    fullfile(project,'script'),fullfile(project,'Atmospheric data')};
for k=1:numel(P),if isfolder(P{k}),addpath(P{k});end,end

wave=(400:5:2500)';M=nk_data(wave);                         % n-i*k
ID.SiO2=2;ID.ITO=4;ID.DAST=10;ID.Au=11;ID.Ag=13;
assert(size(M,2)>=13,'nk_data must contain material IDs through Ag = 13.');
assert(all(isfinite(M(:,[ID.SiO2 ID.ITO ID.DAST ID.Au ID.Ag])),'all'),...
    'Material data contain NaN/Inf over 400-2500 nm.');

%% Baseline
dAu0=5;dAg0=8;dDAST=100;dITO=40;dSi10=80;dSi20=90;        % nm
n0=2.2;dn_du=3.41;Vscan=(-10:1:10)';

%% Robustness ranges
AuScan=4:1:7;
AgScan=6:2:12;
Si1Scan=70:10:100;
Si2Scan=80:10:110;

%% Spectral weights
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
assert(trapz(wave(iLum),Veye(iLum))>0&&trapz(wave(iSol),S(iSol))>0,'Invalid weighting data.');

%% Study A: Au/Ag robustness
NA=numel(AuScan);NG=numel(AgScan);
VcMetal=nan(NG,NA);MaxDAMetal=zeros(NG,NA);Tlum0Metal=zeros(NG,NA);
DTsolMetal=zeros(NG,NA);SpanMetal=zeros(NG,NA);RecMetal=zeros(NG,NA);

fprintf('\nScanning Au/Ag robustness...\n');
for ia=1:NA
    for ig=1:NG
        [VcMetal(ig,ia),MaxDAMetal(ig,ia),Tlum0Metal(ig,ia),DTsolMetal(ig,ia),...
            SpanMetal(ig,ia),RecMetal(ig,ia)]=screenDesign(wave,M,Vscan,...
            AuScan(ia),AgScan(ig),dDAST,dITO,dSi10,dSi20,ID,n0,dn_du,S,Veye,iSol,iLum);
    end
    fprintf('Au = %g nm completed\n',AuScan(ia));
end

%% Study B: SiO2 spacer robustness
N1=numel(Si1Scan);N2=numel(Si2Scan);
VcSi=nan(N2,N1);MaxDASi=zeros(N2,N1);Tlum0Si=zeros(N2,N1);
DTsolSi=zeros(N2,N1);SpanSi=zeros(N2,N1);RecSi=zeros(N2,N1);

fprintf('\nScanning SiO2 spacer robustness...\n');
for i1=1:N1
    for i2=1:N2
        [VcSi(i2,i1),MaxDASi(i2,i1),Tlum0Si(i2,i1),DTsolSi(i2,i1),...
            SpanSi(i2,i1),RecSi(i2,i1)]=screenDesign(wave,M,Vscan,...
            dAu0,dAg0,dDAST,dITO,Si1Scan(i1),Si2Scan(i2),ID,n0,dn_du,S,Veye,iSol,iLum);
    end
    fprintf('SiO2_1 = %g nm completed\n',Si1Scan(i1));
end

%% Baseline indices
[~,iAu0]=min(abs(AuScan-dAu0));[~,iAg0]=min(abs(AgScan-dAg0));
[~,iSi10]=min(abs(Si1Scan-dSi10));[~,iSi20]=min(abs(Si2Scan-dSi20));

%% Summary
fprintf('\n================ ROBUSTNESS SUMMARY =================\n');
fprintf('Baseline: Au=%g nm, Ag=%g nm, SiO2_1=%g nm, SiO2_2=%g nm\n',...
    dAu0,dAg0,dSi10,dSi20);

fprintf('\nMetal-thickness robustness:\n');
fprintf('  Reversal designs       = %d / %d\n',sum(isfinite(VcMetal),'all'),numel(VcMetal));
fprintf('  Vc range               = %.2f to %.2f V\n',min(VcMetal,[],'all','omitnan'),max(VcMetal,[],'all','omitnan'));
fprintf('  Tlum(0V) range         = %.4f to %.4f\n',min(Tlum0Metal,[],'all'),max(Tlum0Metal,[],'all'));
fprintf('  DeltaTsol range        = %.4f to %.4f\n',min(DTsolMetal,[],'all'),max(DTsolMetal,[],'all'));
fprintf('  max|DeltaAsol| range   = %.4f to %.4f\n',min(MaxDAMetal,[],'all'),max(MaxDAMetal,[],'all'));
fprintf('  max reciprocity error  = %.3e\n',max(RecMetal,[],'all'));

fprintf('\nSpacer robustness:\n');
fprintf('  Reversal designs       = %d / %d\n',sum(isfinite(VcSi),'all'),numel(VcSi));
fprintf('  Vc range               = %.2f to %.2f V\n',min(VcSi,[],'all','omitnan'),max(VcSi,[],'all','omitnan'));
fprintf('  Tlum(0V) range         = %.4f to %.4f\n',min(Tlum0Si,[],'all'),max(Tlum0Si,[],'all'));
fprintf('  DeltaTsol range        = %.4f to %.4f\n',min(DTsolSi,[],'all'),max(DTsolSi,[],'all'));
fprintf('  max|DeltaAsol| range   = %.4f to %.4f\n',min(MaxDASi,[],'all'),max(MaxDASi,[],'all'));
fprintf('  max reciprocity error  = %.3e\n',max(RecSi,[],'all'));

fprintf('\nBaseline metrics:\n');
fprintf('  Vc                    = %.3f V\n',VcMetal(iAg0,iAu0));
fprintf('  Tlum(0V)              = %.4f\n',Tlum0Metal(iAg0,iAu0));
fprintf('  DeltaTsol             = %.4f\n',DTsolMetal(iAg0,iAu0));
fprintf('  max|DeltaAsol|        = %.4f\n',MaxDAMetal(iAg0,iAu0));
fprintf('  V1 DeltaAsol span     = %.4f\n',SpanMetal(iAg0,iAu0));
fprintf('=====================================================\n');

%% Compact numerical tables
fprintf('\nAu/Ag robustness table:\n');
fprintf(' Au   Ag    Vc(V)   Tlum0   DeltaTsol   max|DAsol|\n');
for ia=1:NA
    for ig=1:NG
        if isfinite(VcMetal(ig,ia)),vc=sprintf('%+.2f',VcMetal(ig,ia));else,vc='  -- ';end
        fprintf('%3g  %3g   %6s   %.4f    %.4f       %.4f\n',...
            AuScan(ia),AgScan(ig),vc,Tlum0Metal(ig,ia),DTsolMetal(ig,ia),MaxDAMetal(ig,ia));
    end
end

fprintf('\nSiO2 robustness table:\n');
fprintf(' Si1  Si2    Vc(V)   Tlum0   DeltaTsol   max|DAsol|\n');
for i1=1:N1
    for i2=1:N2
        if isfinite(VcSi(i2,i1)),vc=sprintf('%+.2f',VcSi(i2,i1));else,vc='  -- ';end
        fprintf('%3g  %3g   %6s   %.4f    %.4f       %.4f\n',...
            Si1Scan(i1),Si2Scan(i2),vc,Tlum0Si(i2,i1),DTsolSi(i2,i1),MaxDASi(i2,i1));
    end
end

%% Figure
figure('Color','w','Position',[30 50 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');ax=gobjects(1,8);

% (a) Metal Vc
ax(1)=nexttile;
imagesc(ax(1),AuScan,AgScan,VcMetal);axis(ax(1),'xy');caxis(ax(1),[-10 10]);
xlabel(ax(1),'Au thickness (nm)');ylabel(ax(1),'Ag thickness (nm)');
cb=colorbar(ax(1));cb.Label.String='Reversal voltage {\it V}_c (V)';cb.Label.Interpreter='tex';
hold(ax(1),'on');plot(ax(1),dAu0,dAg0,'kp','MarkerSize',11,'LineWidth',1.6);

% (b) Metal directional contrast
ax(2)=nexttile;
imagesc(ax(2),AuScan,AgScan,MaxDAMetal);axis(ax(2),'xy');
xlabel(ax(2),'Au thickness (nm)');ylabel(ax(2),'Ag thickness (nm)');
cb=colorbar(ax(2));cb.Label.String='max|\Delta {\it A}_{sol}|';cb.Label.Interpreter='tex';
hold(ax(2),'on');plot(ax(2),dAu0,dAg0,'kp','MarkerSize',11,'LineWidth',1.6);

% (c) Metal Tlum0
ax(3)=nexttile;
imagesc(ax(3),AuScan,AgScan,Tlum0Metal);axis(ax(3),'xy');caxis(ax(3),[0 1]);
xlabel(ax(3),'Au thickness (nm)');ylabel(ax(3),'Ag thickness (nm)');
cb=colorbar(ax(3));cb.Label.String='{\it T}_{lum}(0 V)';cb.Label.Interpreter='tex';
hold(ax(3),'on');plot(ax(3),dAu0,dAg0,'kp','MarkerSize',11,'LineWidth',1.6);

% (d) Metal DeltaTsol
ax(4)=nexttile;
imagesc(ax(4),AuScan,AgScan,DTsolMetal);axis(ax(4),'xy');
xlabel(ax(4),'Au thickness (nm)');ylabel(ax(4),'Ag thickness (nm)');
cb=colorbar(ax(4));cb.Label.String='\Delta {\it T}_{sol}';cb.Label.Interpreter='tex';
hold(ax(4),'on');plot(ax(4),dAu0,dAg0,'kp','MarkerSize',11,'LineWidth',1.6);

% (e) Spacer Vc
ax(5)=nexttile;
imagesc(ax(5),Si1Scan,Si2Scan,VcSi);axis(ax(5),'xy');caxis(ax(5),[-10 10]);
xlabel(ax(5),'SiO_2 spacer 1 (nm)');ylabel(ax(5),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(5));cb.Label.String='Reversal voltage {\it V}_c (V)';cb.Label.Interpreter='tex';
hold(ax(5),'on');plot(ax(5),dSi10,dSi20,'kp','MarkerSize',11,'LineWidth',1.6);

% (f) Spacer directional contrast
ax(6)=nexttile;
imagesc(ax(6),Si1Scan,Si2Scan,MaxDASi);axis(ax(6),'xy');
xlabel(ax(6),'SiO_2 spacer 1 (nm)');ylabel(ax(6),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(6));cb.Label.String='max|\Delta {\it A}_{sol}|';cb.Label.Interpreter='tex';
hold(ax(6),'on');plot(ax(6),dSi10,dSi20,'kp','MarkerSize',11,'LineWidth',1.6);

% (g) Spacer Tlum0
ax(7)=nexttile;
imagesc(ax(7),Si1Scan,Si2Scan,Tlum0Si);axis(ax(7),'xy');caxis(ax(7),[0 1]);
xlabel(ax(7),'SiO_2 spacer 1 (nm)');ylabel(ax(7),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(7));cb.Label.String='{\it T}_{lum}(0 V)';cb.Label.Interpreter='tex';
hold(ax(7),'on');plot(ax(7),dSi10,dSi20,'kp','MarkerSize',11,'LineWidth',1.6);

% (h) Spacer DeltaTsol
ax(8)=nexttile;
imagesc(ax(8),Si1Scan,Si2Scan,DTsolSi);axis(ax(8),'xy');
xlabel(ax(8),'SiO_2 spacer 1 (nm)');ylabel(ax(8),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(8));cb.Label.String='\Delta {\it T}_{sol}';cb.Label.Interpreter='tex';
hold(ax(8),'on');plot(ax(8),dSi10,dSi20,'kp','MarkerSize',11,'LineWidth',1.6);

Panel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
for j=1:8
    set(ax(j),'FontName',FN,'FontSize',FS1,'FontWeight','normal','LineWidth',LW,...
        'TickDir','in','TickLength',[.012 .012],'Box','on','Layer','top');
    text(ax(j),0,1.06,Panel{j},'Units','normalized','FontName',FN,'FontSize',FS2,...
        'FontWeight','normal','HorizontalAlignment','left','VerticalAlignment','bottom',...
        'Clipping','off','HandleVisibility','off');
end

fprintf('Elapsed time: %.2f s\n',toc);

%% Local functions
function [Vc,MaxDA,Tlum0,DTsol,Span,Rec]=screenDesign(wave,M,Vscan,dAu,dAg,dDAST,dITO,dSi1,dSi2,ID,n0,dn_du,S,Veye,iSol,iLum)
Nv=numel(Vscan);DA=zeros(Nv,1);Tsol=zeros(Nv,1);Rec=zeros(Nv,1);Tlum0=NaN;
d_nm=[dAu dDAST dITO dSi1 dITO dDAST dAg dSi2 dAg dDAST dITO];

for iv=1:Nv
    V=Vscan(iv);

    % V1-only directional response
    [~,Tf1,Af1,~,Tb1,Ab1]=calcState(wave,M,[V 0 0],d_nm,dDAST,ID,n0,dn_du);
    DA(iv)=weighted(Af1-Ab1,wave,S,iSol);

    % Common-voltage smart-window response
    [~,Tf,~,~,Tb,~]=calcState(wave,M,[V V V],d_nm,dDAST,ID,n0,dn_du);
    Tsol(iv)=weighted(Tf,wave,S,iSol);
    Rec(iv)=max(abs(Tf-Tb));

    if V==0
        Tlum0=weighted(Tf,wave,Veye,iLum);
    end
end

Vc=zeroCross(Vscan,DA);
MaxDA=max(abs(DA));
Span=max(DA)-min(DA);
DTsol=max(Tsol)-min(Tsol);
Rec=max(Rec);
end

function [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,V,d_nm,tDAST,ID,n0,dn_du)
N=numel(wave);Rf=zeros(N,1);Tf=Rf;Af=Rf;Rb=Rf;Tb=Rf;Ab=Rf;
for iw=1:N
    d1=dastIndex(M(iw,ID.DAST),V(1),tDAST,n0,dn_du);
    d2=dastIndex(M(iw,ID.DAST),V(2),tDAST,n0,dn_du);
    d3=dastIndex(M(iw,ID.DAST),V(3),tDAST,n0,dn_du);

    nf=[M(iw,ID.Au),d1,M(iw,ID.ITO),M(iw,ID.SiO2),M(iw,ID.ITO),...
        d2,M(iw,ID.Ag),M(iw,ID.SiO2),M(iw,ID.Ag),d3,M(iw,ID.ITO)];

    lam=wave(iw)*1e-9;df=d_nm*1e-9;
    [Rf(iw),Tf(iw)]=tmmRT(1,1,nf,df,lam);
    [Rb(iw),Tb(iw)]=tmmRT(1,1,fliplr(nf),fliplr(df),lam);

    Af(iw)=1-Rf(iw)-Tf(iw);
    Ab(iw)=1-Rb(iw)-Tb(iw);
end
end

function nV=dastIndex(nk0,V,t_nm,n0,dn_du)
nV=(n0+dn_du*V/t_nm)-1i*abs(imag(nk0));
end

function [R,T]=tmmRT(n0,ns,nf,df,lam)
M=eye(2);
for j=1:numel(nf)
    n=nf(j);delta=2*pi*n*df(j)/lam;
    M=M*[cos(delta) 1i*sin(delta)/n;1i*n*sin(delta) cos(delta)];
end
den=n0*M(1,1)+n0*ns*M(1,2)+M(2,1)+ns*M(2,2);
r=(n0*M(1,1)+n0*ns*M(1,2)-M(2,1)-ns*M(2,2))/den;
t=2*n0/den;
R=real(abs(r)^2);
T=real(max(real(ns),0)/max(real(n0),eps)*abs(t)^2);
end

function y=weighted(x,w,W,idx)
x=x(:);w=w(:);W=W(:);
y=trapz(w(idx),x(idx).*W(idx))/trapz(w(idx),W(idx));
end

function Vc=zeroCross(V,y)
Vc=NaN;
for k=1:numel(V)-1
    if y(k)==0,Vc=V(k);return;end
    if y(k)*y(k+1)<0
        Vc=V(k)-y(k)*(V(k+1)-V(k))/(y(k+1)-y(k));
        return;
    end
end
if y(end)==0,Vc=V(end);end
end

function file=findData(root,names)
project=fileparts(root);
folders={root,project,fullfile(root,'Atmospheric data'),fullfile(project,'Atmospheric data')};
file='';
for j=1:numel(names)
    for k=1:numel(folders)
        f=fullfile(folders{k},names{j});
        if isfile(f),file=f;return;end
    end
    f=which(names{j});
    if ~isempty(f),file=f;return;end
end
error('Cannot find data file: %s',strjoin(names,', '));
end

function [w,y]=readWeight(file,col)
M=readmatrix(file);
assert(size(M,2)>=col,'Insufficient columns in %s.',file);
ok=isfinite(M(:,1))&isfinite(M(:,col))&M(:,1)>0;
w=M(ok,1);y=M(ok,col);
if max(w)<100,w=w*1000;end
[w,id]=sort(w);y=y(id);
[w,id]=unique(w,'stable');y=y(id);
end