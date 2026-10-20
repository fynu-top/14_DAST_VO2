% exportgraphics(gcf,'DAST_Sellmeier_validation.pdf','ContentType','vector');
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
    'Material optical constants contain NaN/Inf.');

%% Fixed nominal structure
% Air/Au/DAST1/ITO/SiO2_1/ITO/DAST2/Ag/SiO2_2/Ag/DAST3/ITO/Air
dAu=5;dDAST=100;dITO=40;dSi1=80;dSi2=90;dAg=8;             % nm
d_nm=[dAu dDAST dITO dSi1 dITO dDAST dAg dSi2 dAg dDAST dITO];

%% Electro-optic and Sellmeier parameters
nConst=2.2;dn_du=3.41;                                      % nm/V
qX=1.645;nInfX=2.078;lambda0X=533;                          % nm
BlendLow=650;BlendHigh=750;                                 % smooth edge treatment
Vscan=(-10:1:10)';scaleScan=(0.8:0.1:1.2)';
ModelName={'Constant {\it n}_0','Sellmeier-hybrid {\it n}_x(\lambda)'};

%% Spectral weights
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
DenSol=trapz(wave(iSol),S(iSol));
assert(DenSol>0&&trapz(wave(iLum),Veye(iLum))>0,'Invalid spectral weights.');

%% Show the actual DAST baseline models
nCurrent=nConst*ones(size(wave));
nHybrid=dastBaseN(wave,nConst,qX,nInfX,lambda0X,BlendLow,BlendHigh);

fprintf('\n================ DAST SELMEIER MODEL =================\n');
fprintf('Sellmeier x-axis: q=%.3f, n_inf=%.3f, lambda0=%.1f nm\n',qX,nInfX,lambda0X);
fprintf('Hybrid transition: %.0f-%.0f nm\n',BlendLow,BlendHigh);
for lam=[700 800 1313 1550 1900 2500]
    nRaw=sellmeierX(lam,qX,nInfX,lambda0X);
    nHy=interp1(wave,nHybrid,lam,'linear');
    fprintf('lambda=%4g nm: raw Sellmeier n=%.4f, hybrid n=%.4f\n',lam,nRaw,nHy);
end
fprintf('Current DAST k range = %.4g - %.4g\n',...
    min(abs(imag(M(:,ID.DAST)))),max(abs(imag(M(:,ID.DAST)))));
fprintf('======================================================\n');

%% Compare constant and Sellmeier-hybrid models
Nm=2;Nv=numel(Vscan);
Tlum=zeros(Nv,Nm);Tsol=Tlum;DAsol=Tlum;Rec=zeros(1,Nm);Vc=nan(1,Nm);

fprintf('\nCalculating constant-vs-dispersive comparison...\n');
for m=1:Nm
    [Tlum(:,m),Tsol(:,m),DAsol(:,m),Rec(m)]=calcMetrics(wave,M,Vscan,d_nm,dDAST,ID,...
        nConst,dn_du,1,m,qX,nInfX,lambda0X,BlendLow,BlendHigh,S,Veye,iSol,iLum);
    Vc(m)=findVc(Vscan,DAsol(:,m),7.23);
    fprintf('%s completed\n',ModelName{m});
end

%% Refine Sellmeier-hybrid reversal voltage
fun=@(v)solarDirAtV(v,wave,M,d_nm,dDAST,ID,nConst,dn_du,1,2,...
    qX,nInfX,lambda0X,BlendLow,BlendHigh,S,iSol);
VcSell=refineRoot(fun,Vscan,DAsol(:,2),Vc(2));
Vc(2)=VcSell;

%% Mechanism spectra for Sellmeier-hybrid model
Vstate=[-10 VcSell 10];
StateName={sprintf('{\\it V}_1=-10 V'),sprintf('{\\it V}_1={\\it V}_c=%.2f V',VcSell),...
           sprintf('{\\it V}_1=+10 V')};
DAstate=zeros(numel(wave),3);Cspec=zeros(numel(wave),3);Ccum=zeros(numel(wave),3);
DAsolState=zeros(1,3);VisContribution=zeros(1,3);NIRContribution=zeros(1,3);

for s=1:3
    [~,~,Af,~,~,Ab]=calcState(wave,M,[Vstate(s) 0 0],d_nm,dDAST,ID,nConst,dn_du,1,2,...
        qX,nInfX,lambda0X,BlendLow,BlendHigh);
    DAstate(:,s)=Af-Ab;
    Cspec(:,s)=S.*DAstate(:,s)/DenSol;
    Ccum(:,s)=cumtrapz(wave,Cspec(:,s));
    DAsolState(s)=Ccum(end,s);
    ivis=wave>=400&wave<700;inir=wave>=700&wave<=2500;
    VisContribution(s)=trapz(wave(ivis),Cspec(ivis,s));
    NIRContribution(s)=trapz(wave(inir),Cspec(inir,s));
end

%% EO-coefficient sensitivity using Sellmeier-hybrid model
Ns=numel(scaleScan);VcScale=nan(Ns,1);DTlumScale=zeros(Ns,1);
DTsolScale=zeros(Ns,1);MaxDAScale=zeros(Ns,1);RecScale=zeros(Ns,1);

fprintf('Calculating Sellmeier-hybrid EO-coefficient sensitivity...\n');
for is=1:Ns
    sc=scaleScan(is);
    [tl,ts,da,RecScale(is)]=calcMetrics(wave,M,Vscan,d_nm,dDAST,ID,nConst,dn_du,sc,2,...
        qX,nInfX,lambda0X,BlendLow,BlendHigh,S,Veye,iSol,iLum);
    VcScale(is)=findVc(Vscan,da,VcSell);
    DTlumScale(is)=max(tl)-min(tl);
    DTsolScale(is)=max(ts)-min(ts);
    MaxDAScale(is)=max(abs(da));
end

%% Console summary
fprintf('\n================ DISPERSION VALIDATION =================\n');
fprintf('Model                         Vc(V)    DeltaTlum   DeltaTsol   max|DeltaAsol|   max|Tf-Tb|\n');
for m=1:Nm
    fprintf('%-28s %+7.3f     %.4f       %.4f        %.4f        %.3e\n',...
        ModelName{m},Vc(m),max(Tlum(:,m))-min(Tlum(:,m)),...
        max(Tsol(:,m))-min(Tsol(:,m)),max(abs(DAsol(:,m))),Rec(m));
end

fprintf('\nSellmeier-hybrid directional reversal:\n');
for s=1:3
    fprintf('V1=%+7.3f V: DeltaAsol=%+.5f, 400-700 nm contribution=%+.5f, 700-2500 nm contribution=%+.5f\n',...
        Vstate(s),DAsolState(s),VisContribution(s),NIRContribution(s));
end

fprintf('\nSellmeier-hybrid EO sensitivity:\n');
fprintf(' scale  dn/du(nm/V)    Vc(V)   DeltaTlum  DeltaTsol  max|DeltaAsol|\n');
for is=1:Ns
    fprintf(' %.1f      %.3f       %+6.2f    %.4f      %.4f       %.4f\n',...
        scaleScan(is),dn_du*scaleScan(is),VcScale(is),DTlumScale(is),...
        DTsolScale(is),MaxDAScale(is));
end
fprintf('Maximum reciprocity error = %.3e\n',max([Rec(:);RecScale(:)]));
fprintf('========================================================\n');

%% Figure
figure('Color','w','Position',[30 50 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
ax=gobjects(1,8);C=lines(3);

% (a) DAST dispersion
ax(1)=nexttile;hold(ax(1),'on');
plot(ax(1),wave/1000,nCurrent,'--','LineWidth',LW);
plot(ax(1),wave/1000,nHybrid,'-','LineWidth',LW);
xline(ax(1),.7,':','LineWidth',1,'HandleVisibility','off');
xlabel(ax(1),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(1),'{\it n}_{DAST}','Interpreter','tex');
legend(ax(1),ModelName,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(1),[.4 2.5]);

% (b) Common-voltage Tlum
ax(2)=nexttile;hold(ax(2),'on');
for m=1:Nm,plot(ax(2),Vscan,Tlum(:,m),'LineWidth',LW);end
xlabel(ax(2),'Common voltage (V)');
ylabel(ax(2),'{\it T}_{lum}','Interpreter','tex');
legend(ax(2),ModelName,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(2),[-10 10]);ylim(ax(2),[0 1]);

% (c) Common-voltage Tsol
ax(3)=nexttile;hold(ax(3),'on');
for m=1:Nm,plot(ax(3),Vscan,Tsol(:,m),'LineWidth',LW);end
xlabel(ax(3),'Common voltage (V)');
ylabel(ax(3),'{\it T}_{sol}','Interpreter','tex');
legend(ax(3),ModelName,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(3),[-10 10]);ylim(ax(3),[0 1]);

% (d) V1-only reversal
ax(4)=nexttile;hold(ax(4),'on');
for m=1:Nm,plot(ax(4),Vscan,DAsol(:,m),'LineWidth',LW);end
yline(ax(4),0,'k:','HandleVisibility','off');
xline(ax(4),VcSell,':','LineWidth',1,'HandleVisibility','off');
xlabel(ax(4),'{\it V}_1 only (V)','Interpreter','tex');
ylabel(ax(4),'\Delta {\it A}_{sol}','Interpreter','tex');
legend(ax(4),ModelName,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(4),[-10 10]);

% (e) Sellmeier spectral directional absorption
ax(5)=nexttile;hold(ax(5),'on');
for s=1:3,plot(ax(5),wave/1000,DAstate(:,s),'LineWidth',LW,'Color',C(s,:));end
yline(ax(5),0,'k:','HandleVisibility','off');
xlabel(ax(5),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(5),'\Delta {\it A}_{dir}','Interpreter','tex');
legend(ax(5),StateName,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(5),[.4 2.5]);

% (f) Solar-weighted spectral contribution density
ax(6)=nexttile;hold(ax(6),'on');
for s=1:3,plot(ax(6),wave/1000,Cspec(:,s)*1000,'LineWidth',LW,'Color',C(s,:));end
yline(ax(6),0,'k:','HandleVisibility','off');
xlabel(ax(6),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(6),'Solar-weighted contribution (10^{-3} nm^{-1})');
legend(ax(6),StateName,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(6),[.4 2.5]);

% (g) Cumulative directional contribution
ax(7)=nexttile;hold(ax(7),'on');
for s=1:3,plot(ax(7),wave/1000,Ccum(:,s),'LineWidth',LW,'Color',C(s,:));end
yline(ax(7),0,'k:','HandleVisibility','off');
xlabel(ax(7),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(7),'Cumulative \Delta {\it A}_{sol}','Interpreter','tex');
legend(ax(7),StateName,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(7),[.4 2.5]);

% (h) EO coefficient robustness under dispersive model
ax(8)=nexttile;hold(ax(8),'on');
yyaxis(ax(8),'left');
plot(ax(8),scaleScan,VcScale,'-o','LineWidth',LW,'MarkerSize',6);
ylabel(ax(8),'Reversal voltage {\it V}_c (V)','Interpreter','tex');
yyaxis(ax(8),'right');
plot(ax(8),scaleScan,DTsolScale,'-s','LineWidth',LW,'MarkerSize',6);
plot(ax(8),scaleScan,MaxDAScale,'-^','LineWidth',LW,'MarkerSize',6);
ylabel(ax(8),'Optical metric');
xline(ax(8),1,':','HandleVisibility','off');
xlabel(ax(8),'Electro-optic coefficient scale');
legend(ax(8),{'{\it V}_c','\Delta {\it T}_{sol}','max|\Delta {\it A}_{sol}|'},...
    'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(8),[.8 1.2]);

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
function [Tlum,Tsol,DAsol,Rec]=calcMetrics(wave,M,Vscan,d_nm,tDAST,ID,nConst,dn_du,...
    scale,mode,qX,nInfX,lambda0X,b1,b2,S,Veye,iSol,iLum)
Nv=numel(Vscan);Tlum=zeros(Nv,1);Tsol=Tlum;DAsol=Tlum;Rec=0;
for iv=1:Nv
    V=Vscan(iv);
    [~,Tf,~,~,Tb,~]=calcState(wave,M,[V V V],d_nm,tDAST,ID,nConst,dn_du,scale,mode,...
        qX,nInfX,lambda0X,b1,b2);
    Tlum(iv)=weighted(Tf,wave,Veye,iLum);
    Tsol(iv)=weighted(Tf,wave,S,iSol);
    Rec=max(Rec,max(abs(Tf-Tb)));

    [~,~,Af,~,~,Ab]=calcState(wave,M,[V 0 0],d_nm,tDAST,ID,nConst,dn_du,scale,mode,...
        qX,nInfX,lambda0X,b1,b2);
    DAsol(iv)=weighted(Af-Ab,wave,S,iSol);
end
end

function [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,V,d_nm,tDAST,ID,nConst,dn_du,...
    scale,mode,qX,nInfX,lambda0X,b1,b2)
N=numel(wave);Rf=zeros(N,1);Tf=Rf;Af=Rf;Rb=Rf;Tb=Rf;Ab=Rf;
for iw=1:N
    d1=dastIndex(wave(iw),M(iw,ID.DAST),V(1),tDAST,nConst,dn_du,scale,mode,qX,nInfX,lambda0X,b1,b2);
    d2=dastIndex(wave(iw),M(iw,ID.DAST),V(2),tDAST,nConst,dn_du,scale,mode,qX,nInfX,lambda0X,b1,b2);
    d3=dastIndex(wave(iw),M(iw,ID.DAST),V(3),tDAST,nConst,dn_du,scale,mode,qX,nInfX,lambda0X,b1,b2);

    nf=[M(iw,ID.Au),d1,M(iw,ID.ITO),M(iw,ID.SiO2),M(iw,ID.ITO),...
        d2,M(iw,ID.Ag),M(iw,ID.SiO2),M(iw,ID.Ag),d3,M(iw,ID.ITO)];

    lam=wave(iw)*1e-9;df=d_nm*1e-9;
    [Rf(iw),Tf(iw)]=tmmRT(1,1,nf,df,lam);
    [Rb(iw),Tb(iw)]=tmmRT(1,1,fliplr(nf),fliplr(df),lam);
    Af(iw)=1-Rf(iw)-Tf(iw);Ab(iw)=1-Rb(iw)-Tb(iw);
end
end

function nV=dastIndex(w,nk0,V,t_nm,nConst,dn_du,scale,mode,qX,nInfX,lambda0X,b1,b2)
if mode==1
    nr=nConst;
else
    nr=dastBaseN(w,nConst,qX,nInfX,lambda0X,b1,b2);
end
nV=(nr+scale*dn_du*V/t_nm)-1i*abs(imag(nk0));
end

function n=dastBaseN(w,nConst,q,nInf,l0,b1,b2)
n=nConst*ones(size(w));
id2=w>=b2;
n(id2)=sellmeierX(w(id2),q,nInf,l0);
id=w>b1&w<b2;
if any(id)
    s=.5-.5*cos(pi*(w(id)-b1)/(b2-b1));
    ns=sellmeierX(w(id),q,nInf,l0);
    n(id)=(1-s).*nConst+s.*ns;
end
end

function n=sellmeierX(w,q,nInf,l0)
n=sqrt(nInf.^2+q*l0.^2./(w.^2-l0.^2));
end

function val=solarDirAtV(V1,wave,M,d_nm,tDAST,ID,nConst,dn_du,scale,mode,...
    qX,nInfX,l0,b1,b2,S,iSol)
[~,~,Af,~,~,Ab]=calcState(wave,M,[V1 0 0],d_nm,tDAST,ID,nConst,dn_du,...
    scale,mode,qX,nInfX,l0,b1,b2);
val=weighted(Af-Ab,wave,S,iSol);
end

function Vc=findVc(V,y,Vref)
roots=[];
for k=1:numel(V)-1
    if y(k)==0,roots(end+1)=V(k);end
    if y(k)*y(k+1)<0
        roots(end+1)=V(k)-y(k)*(V(k+1)-V(k))/(y(k+1)-y(k));
    end
end
if y(end)==0,roots(end+1)=V(end);end
if isempty(roots),Vc=NaN;else,[~,i]=min(abs(roots-Vref));Vc=roots(i);end
end

function Vc=refineRoot(fun,V,y,Vguess)
Vc=Vguess;
for k=1:numel(V)-1
    if y(k)*y(k+1)<0
        r=fzero(fun,[V(k) V(k+1)]);
        if isnan(Vc)||abs(r-Vguess)<abs(Vc-Vguess),Vc=r;end
    end
end
end

function [R,T]=tmmRT(n0,ns,nf,df,lam)
M=eye(2);
for j=1:numel(nf)
    n=nf(j);d=2*pi*n*df(j)/lam;
    M=M*[cos(d) 1i*sin(d)/n;1i*n*sin(d) cos(d)];
end
den=n0*M(1,1)+n0*ns*M(1,2)+M(2,1)+ns*M(2,2);
r=(n0*M(1,1)+n0*ns*M(1,2)-M(2,1)-ns*M(2,2))/den;t=2*n0/den;
R=real(abs(r)^2);T=real(max(real(ns),0)/max(real(n0),eps)*abs(t)^2);
end

function y=weighted(x,w,W,idx)
y=trapz(w(idx),x(idx).*W(idx))/trapz(w(idx),W(idx));
end

function file=findData(root,names)
project=fileparts(root);
folders={root,project,fullfile(root,'Atmospheric data'),fullfile(project,'Atmospheric data')};
file='';
for j=1:numel(names)
    for k=1:numel(folders)
        f=fullfile(folders{k},names{j});if isfile(f),file=f;return;end
    end
    f=which(names{j});if ~isempty(f),file=f;return;end
end
error('Cannot find data file: %s',strjoin(names,', '));
end

function [w,y]=readWeight(file,col)
M=readmatrix(file);ok=isfinite(M(:,1))&isfinite(M(:,col))&M(:,1)>0;
w=M(ok,1);y=M(ok,col);if max(w)<100,w=w*1000;end
[w,id]=sort(w);y=y(id);[w,id]=unique(w,'stable');y=y(id);
end