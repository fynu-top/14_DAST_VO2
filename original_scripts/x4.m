% exportgraphics(gcf,'DAST_core_performance.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

FS1=20;FS2=26;FN='Times New Roman';LW=1.5;FigWidth=1900;FigHeight=880;mu=char(181);
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

wave=(400:5:2500)';M=nk_data(wave);                         % n-i*k convention
ID.SiO2=2;ID.ITO=4;ID.DAST=10;ID.Au=11;ID.Ag=13;
assert(size(M,2)>=13,'nk_data must contain material IDs through Ag = 13.');
assert(all(isfinite(M(:,[ID.SiO2 ID.ITO ID.DAST ID.Au ID.Ag])),'all'),...
    'Material optical constants contain NaN/Inf over 400-2500 nm.');

%% Final baseline structure
% Air/Au/DAST1/ITO/SiO2_1/ITO/DAST2/Ag/SiO2_2/Ag/DAST3/ITO/Air
dAu=5;dDAST=100;dITO=40;dSi1=80;dSi2=90;dAg=8;             % nm
d_nm=[dAu dDAST dITO dSi1 dITO dDAST dAg dSi2 dAg dDAST dITO];

n0=2.2;dn_du=3.41;                                          % nm/V
Vscan=(-10:1:10)';Vshow=[-10 0 10];
Nv=numel(Vscan);Nw=numel(wave);

fprintf('\n================ FINAL STRUCTURE =================\n');
fprintf('Air/Au(5)/DAST(100)/ITO(40)/SiO2(80)/ITO(40)/DAST(100)/Ag(8)/SiO2(90)/Ag(8)/DAST(100)/ITO(40)/Air\n');
fprintf('DAST: n(-10V)=%.4f, n(0V)=%.4f, n(+10V)=%.4f\n',...
    n0-dn_du*10/dDAST,n0,n0+dn_du*10/dDAST);
fprintf('|Delta n|max = %.4f, |E|max = %.1f V/um\n',dn_du*10/dDAST,10000/dDAST);
fprintf('==================================================\n');

%% Spectral weighting
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
assert(trapz(wave(iLum),Veye(iLum))>0&&trapz(wave(iSol),S(iSol))>0,'Invalid spectral weights.');

%% Common-voltage response
RfMap=zeros(Nv,Nw);TfMap=RfMap;AfMap=RfMap;RbMap=RfMap;TbMap=RfMap;AbMap=RfMap;
TlumC=zeros(Nv,1);TsolC=TlumC;AsolFC=TlumC;AsolBC=TlumC;DAsolC=TlumC;
RecErr=TlumC;EnergyErr=TlumC;

fprintf('\nCalculating common-voltage response...\n');
for iv=1:Nv
    V=Vscan(iv);
    [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,[V V V],d_nm,dDAST,ID,n0,dn_du);

    RfMap(iv,:)=Rf;TfMap(iv,:)=Tf;AfMap(iv,:)=Af;
    RbMap(iv,:)=Rb;TbMap(iv,:)=Tb;AbMap(iv,:)=Ab;

    TlumC(iv)=weighted(Tf,wave,Veye,iLum);
    TsolC(iv)=weighted(Tf,wave,S,iSol);
    AsolFC(iv)=weighted(Af,wave,S,iSol);
    AsolBC(iv)=weighted(Ab,wave,S,iSol);
    DAsolC(iv)=AsolFC(iv)-AsolBC(iv);

    RecErr(iv)=max(abs(Tf-Tb));
    EnergyErr(iv)=max([abs(Rf+Tf+Af-1);abs(Rb+Tb+Ab-1)]);
end

%% Independent V1/V2/V3 control
TlumInd=zeros(Nv,3);TsolInd=TlumInd;DAsolInd=TlumInd;
VcInd=nan(1,3);DTsolInd=zeros(1,3);DTlumInd=zeros(1,3);MaxDAInd=zeros(1,3);

fprintf('Calculating independent V1/V2/V3 response...\n');
for j=1:3
    for iv=1:Nv
        VV=[0 0 0];VV(j)=Vscan(iv);
        [~,Tf,Af,~,~,Ab]=calcState(wave,M,VV,d_nm,dDAST,ID,n0,dn_du);
        TlumInd(iv,j)=weighted(Tf,wave,Veye,iLum);
        TsolInd(iv,j)=weighted(Tf,wave,S,iSol);
        DAsolInd(iv,j)=weighted(Af,wave,S,iSol)-weighted(Ab,wave,S,iSol);
    end
    VcInd(j)=zeroCross(Vscan,DAsolInd(:,j));
    DTsolInd(j)=max(TsolInd(:,j))-min(TsolInd(:,j));
    DTlumInd(j)=max(TlumInd(:,j))-min(TlumInd(:,j));
    MaxDAInd(j)=max(abs(DAsolInd(:,j)));
end

VcCommon=zeroCross(Vscan,DAsolC);
DTsolCommon=max(TsolC)-min(TsolC);
DTlumCommon=max(TlumC)-min(TlumC);
MaxDACommon=max(abs(DAsolC));

%% Selected-voltage spectra
Sel=struct;
for q=1:3
    iv=find(Vscan==Vshow(q),1);
    Sel(q).Rf=RfMap(iv,:);Sel(q).Tf=TfMap(iv,:);Sel(q).Af=AfMap(iv,:);
    Sel(q).Rb=RbMap(iv,:);Sel(q).Tb=TbMap(iv,:);Sel(q).Ab=AbMap(iv,:);
end

%% Summary
fprintf('\n================ CORE PERFORMANCE =================\n');
fprintf('Common voltage:\n');
fprintf('  Tlum max/min       = %.4f / %.4f\n',max(TlumC),min(TlumC));
fprintf('  Delta Tlum         = %.4f (%.2f percentage points)\n',DTlumCommon,100*DTlumCommon);
fprintf('  Tsol max/min       = %.4f / %.4f\n',max(TsolC),min(TsolC));
fprintf('  Delta Tsol         = %.4f (%.2f percentage points)\n',DTsolCommon,100*DTsolCommon);
fprintf('  max|Delta Asol|    = %.4f\n',MaxDACommon);
if isfinite(VcCommon),fprintf('  directional Vc     = %.2f V\n',VcCommon);
else,fprintf('  directional Vc     = not found in [-10,+10] V\n');end
fprintf('  max|Tf-Tb|         = %.3e\n',max(RecErr));
fprintf('  max energy error   = %.3e\n',max(EnergyErr));

for Vq=Vshow
    iv=find(Vscan==Vq,1);
    fprintf('\nV=%+3g V: Tlum=%.4f, Tsol=%.4f, Asol,F=%.4f, Asol,B=%.4f, DeltaAsol=%+.4f\n',...
        Vq,TlumC(iv),TsolC(iv),AsolFC(iv),AsolBC(iv),DAsolC(iv));
end

fprintf('\nIndependent channels:\n');
for j=1:3
    fprintf('  V%d: DeltaTlum=%.4f, DeltaTsol=%.4f, max|DeltaAsol|=%.4f, ',...
        j,DTlumInd(j),DTsolInd(j),MaxDAInd(j));
    if isfinite(VcInd(j)),fprintf('Vc=%+.2f V\n',VcInd(j));
    else,fprintf('Vc=none\n');end
end
fprintf('===================================================\n');

%% Figure
figure('Color','w','Position',[30 60 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
ax=gobjects(1,8);C=lines(3);

% (a) Common-voltage transmission spectra
ax(1)=nexttile;hold(ax(1),'on');
for q=1:3
    plot(ax(1),wave/1000,Sel(q).Tf,'LineWidth',LW,'Color',C(q,:));
end
xlabel(ax(1),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(1),'{\it T}','Interpreter','tex');
legend(ax(1),{'-10 V','0 V','+10 V'},'Location','best','Box','off');
xlim(ax(1),[.4 2.5]);ylim(ax(1),[0 1]);

% (b) Forward/reverse R/T/A at 0 V
ax(2)=nexttile;hold(ax(2),'on');q=2;
plot(ax(2),wave/1000,Sel(q).Tf,'k-','LineWidth',LW+0.2);
plot(ax(2),wave/1000,Sel(q).Rf,'-','LineWidth',LW);
plot(ax(2),wave/1000,Sel(q).Rb,'--','LineWidth',LW);
plot(ax(2),wave/1000,Sel(q).Af,'-','LineWidth',LW);
plot(ax(2),wave/1000,Sel(q).Ab,'--','LineWidth',LW);
xlabel(ax(2),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(2),'Optical response');
legend(ax(2),{'{\it T}_F={\it T}_B','{\it R}_F','{\it R}_B','{\it A}_F','{\it A}_B'},...
    'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(2),[.4 2.5]);ylim(ax(2),[0 1]);

% (c) T(lambda,V)
ax(3)=nexttile;
imagesc(ax(3),wave/1000,Vscan,TfMap);axis(ax(3),'xy');caxis(ax(3),[0 1]);
xlabel(ax(3),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(3),'Common voltage (V)');
cb=colorbar(ax(3));cb.Label.String='{\it T}';cb.Label.Interpreter='tex';

% (d) Delta Adir(lambda,V)
ax(4)=nexttile;
Adir=AfMap-AbMap;m=max(abs(Adir),[],'all');if m<1e-6,m=1e-6;end
imagesc(ax(4),wave/1000,Vscan,Adir);axis(ax(4),'xy');caxis(ax(4),[-m m]);
xlabel(ax(4),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(4),'Common voltage (V)');
cb=colorbar(ax(4));cb.Label.String='\Delta {\it A}_{dir}={\it A}_F-{\it A}_B';cb.Label.Interpreter='tex';

% (e) Weighted common-voltage response
ax(5)=nexttile;hold(ax(5),'on');
plot(ax(5),Vscan,TlumC,'-o','LineWidth',LW,'MarkerSize',5);
plot(ax(5),Vscan,TsolC,'-s','LineWidth',LW,'MarkerSize',5);
xlabel(ax(5),'Common voltage (V)');
ylabel(ax(5),'Weighted transmittance');
legend(ax(5),{'{\it T}_{lum}','{\it T}_{sol}'},'Location','best','Box','off','Interpreter','tex');
xlim(ax(5),[-10 10]);ylim(ax(5),[0 1]);

% (f) Independent Tsol channels
ax(6)=nexttile;hold(ax(6),'on');
plot(ax(6),Vscan,TsolInd(:,1),'-','LineWidth',LW);
plot(ax(6),Vscan,TsolInd(:,2),'--','LineWidth',LW);
plot(ax(6),Vscan,TsolInd(:,3),'-.','LineWidth',LW);
plot(ax(6),Vscan,TsolC,':','LineWidth',LW+0.4);
xlabel(ax(6),'Applied voltage (V)');
ylabel(ax(6),'{\it T}_{sol}','Interpreter','tex');
legend(ax(6),{'{\it V}_1 only','{\it V}_2 only','{\it V}_3 only','Common voltage'},...
    'Location','best','Box','off','Interpreter','tex');
xlim(ax(6),[-10 10]);

% (g) Independent directional solar response
ax(7)=nexttile;hold(ax(7),'on');
plot(ax(7),Vscan,DAsolInd(:,1),'-','LineWidth',LW);
plot(ax(7),Vscan,DAsolInd(:,2),'--','LineWidth',LW);
plot(ax(7),Vscan,DAsolInd(:,3),'-.','LineWidth',LW);
plot(ax(7),Vscan,DAsolC,':','LineWidth',LW+0.4);
yline(ax(7),0,'k:','LineWidth',1,'HandleVisibility','off');
xlabel(ax(7),'Applied voltage (V)');
ylabel(ax(7),'\Delta {\it A}_{sol}','Interpreter','tex');
legend(ax(7),{'{\it V}_1 only','{\it V}_2 only','{\it V}_3 only','Common voltage'},...
    'Location','best','Box','off','Interpreter','tex');
xlim(ax(7),[-10 10]);

% (h) Directional modulation and reversal voltage
ax(8)=nexttile;hold(ax(8),'on');
X=1:4;
MaxDA=[MaxDAInd MaxDACommon];
VcAll=[VcInd VcCommon];

yyaxis(ax(8),'left');
b=bar(ax(8),X,MaxDA,0.55);
ylabel(ax(8),'max|\Delta {\it A}_{sol}|','Interpreter','tex');

yyaxis(ax(8),'right');
valid=isfinite(VcAll);
plot(ax(8),X(valid),VcAll(valid),'ko','MarkerFaceColor','w','MarkerSize',8,'LineWidth',LW);
yline(ax(8),0,'k:','HandleVisibility','off');
ylabel(ax(8),'Reversal voltage {\it V}_c (V)','Interpreter','tex');
ylim(ax(8),[-10 10]);
set(ax(8),'XTick',X,'XTickLabel',{'{\it V}_1','{\it V}_2','{\it V}_3','Common'},...
    'TickLabelInterpreter','tex');
xlabel(ax(8),'Voltage channel');

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