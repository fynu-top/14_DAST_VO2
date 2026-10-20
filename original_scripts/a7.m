% exportgraphics(gcf,'DAST_model_sensitivity.pdf','ContentType','vector');
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

%% Fixed structure
dAu=5;dDAST=100;dITO=40;dSi1=80;dSi2=90;dAg=8;             % nm
d_nm=[dAu dDAST dITO dSi1 dITO dDAST dAg dSi2 dAg dDAST dITO];
n0=2.2;dn_du=3.41;                                          % nm/V
Vscan=(-10:1:10)';scaleScan=0.8:0.1:1.2;
ModelName={'Constant n_0','Database dispersion','Anchored dispersion'};
Nm=3;Ns=numel(scaleScan);

%% Database dispersion anchor at 1.55 um
[~,iRef]=min(abs(wave-1550));
nDB=real(M(:,ID.DAST));kDB=abs(imag(M(:,ID.DAST)));
anchorShift=n0-nDB(iRef);

fprintf('\n================ DAST OPTICAL MODELS =================\n');
fprintf('Database n range over 0.4-2.5 um = %.4f - %.4f\n',min(nDB),max(nDB));
fprintf('Database k range over 0.4-2.5 um = %.4f - %.4f\n',min(kDB),max(kDB));
fprintf('n_database(1.55 um)              = %.4f\n',nDB(iRef));
fprintf('Anchoring shift                  = %+.4f\n',anchorShift);
fprintf('Current |Delta n| at +/-10 V     = %.4f\n',dn_du*10/dDAST);
fprintf('======================================================\n');

%% Weight functions
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
assert(trapz(wave(iLum),Veye(iLum))>0&&trapz(wave(iSol),S(iSol))>0,'Invalid spectral weights.');

%% Zero-bias optical models for plotting
nPlot=zeros(numel(wave),Nm);
nPlot(:,1)=n0;
nPlot(:,2)=nDB;
nPlot(:,3)=nDB+anchorShift;

%% Nominal electro-optic coefficient: full comparison
Tlum=zeros(numel(Vscan),Nm);Tsol=Tlum;DAsol=Tlum;
Vc=nan(1,Nm);DTlum=zeros(1,Nm);DTsol=zeros(1,Nm);MaxDA=zeros(1,Nm);RecErr=zeros(1,Nm);

fprintf('\nCalculating nominal model comparison...\n');
for m=1:Nm
    [Tlum(:,m),Tsol(:,m),DAsol(:,m),RecErr(m)]=calcMetrics(...
        wave,M,Vscan,d_nm,dDAST,ID,n0,dn_du,1,m,anchorShift,S,Veye,iSol,iLum);
    Vc(m)=zeroCrossNearest(Vscan,DAsol(:,m),7.23);
    DTlum(m)=max(Tlum(:,m))-min(Tlum(:,m));
    DTsol(m)=max(Tsol(:,m))-min(Tsol(:,m));
    MaxDA(m)=max(abs(DAsol(:,m)));
    fprintf('%s completed\n',ModelName{m});
end

%% Sensitivity to electro-optic coefficient
VcScale=nan(Ns,Nm);DTlumScale=zeros(Ns,Nm);DTsolScale=zeros(Ns,Nm);
MaxDAScale=zeros(Ns,Nm);RecScale=zeros(Ns,Nm);

fprintf('\nCalculating electro-optic-coefficient sensitivity...\n');
for m=1:Nm
    refVc=Vc(m);if ~isfinite(refVc),refVc=7.23;end
    for is=1:Ns
        sc=scaleScan(is);
        [tl,ts,da,rr]=calcMetrics(wave,M,Vscan,d_nm,dDAST,ID,n0,dn_du,sc,m,...
            anchorShift,S,Veye,iSol,iLum);
        VcScale(is,m)=zeroCrossNearest(Vscan,da,refVc);
        DTlumScale(is,m)=max(tl)-min(tl);
        DTsolScale(is,m)=max(ts)-min(ts);
        MaxDAScale(is,m)=max(abs(da));
        RecScale(is,m)=rr;
    end
end

%% Summary
fprintf('\n================ MODEL SENSITIVITY SUMMARY =================\n');
fprintf('Model                         Vc(V)   DeltaTlum   DeltaTsol   max|DeltaAsol|   max|Tf-Tb|\n');
for m=1:Nm
    if isfinite(Vc(m)),vc=sprintf('%+.2f',Vc(m));else,vc='  -- ';end
    fprintf('%-28s %6s     %.4f       %.4f        %.4f        %.3e\n',...
        ModelName{m},vc,DTlum(m),DTsol(m),MaxDA(m),RecErr(m));
end

fprintf('\nEO-coefficient sensitivity:\n');
for m=1:Nm
    fprintf('\n%s\n',ModelName{m});
    fprintf(' scale   dn/du(nm/V)    Vc(V)   DeltaTlum   DeltaTsol   max|DeltaAsol|\n');
    for is=1:Ns
        if isfinite(VcScale(is,m)),vc=sprintf('%+.2f',VcScale(is,m));else,vc='  -- ';end
        fprintf(' %.1f       %.3f       %6s     %.4f       %.4f        %.4f\n',...
            scaleScan(is),scaleScan(is)*dn_du,vc,DTlumScale(is,m),...
            DTsolScale(is,m),MaxDAScale(is,m));
    end
end
fprintf('\nMaximum reciprocity error over all tests = %.3e\n',max(RecScale,[],'all'));
fprintf('============================================================\n');

%% Figure
figure('Color','w','Position',[30 50 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
ax=gobjects(1,8);C=lines(Nm);

% (a) DAST zero-bias optical constants
ax(1)=nexttile;hold(ax(1),'on');
yyaxis(ax(1),'left');
for m=1:Nm,plot(ax(1),wave/1000,nPlot(:,m),'LineWidth',LW,'Color',C(m,:));end
ylabel(ax(1),'{\it n}_{DAST}','Interpreter','tex');
yyaxis(ax(1),'right');
plot(ax(1),wave/1000,kDB,'k--','LineWidth',LW);
ylabel(ax(1),'{\it k}_{DAST}','Interpreter','tex');
xlabel(ax(1),['Wavelength (' mu 'm)'],'Interpreter','none');
legend(ax(1),[ModelName {'{\it k}_{DAST}'}],'Location','best','Box','off',...
    'Interpreter','tex','FontSize',FS1-4);
xlim(ax(1),[.4 2.5]);

% (b) Tlum(V)
ax(2)=nexttile;hold(ax(2),'on');
for m=1:Nm,plot(ax(2),Vscan,Tlum(:,m),'LineWidth',LW,'Color',C(m,:));end
xlabel(ax(2),'Common voltage (V)');
ylabel(ax(2),'{\it T}_{lum}','Interpreter','tex');
legend(ax(2),ModelName,'Location','best','Box','off','FontSize',FS1-4);
xlim(ax(2),[-10 10]);ylim(ax(2),[0 1]);

% (c) Tsol(V)
ax(3)=nexttile;hold(ax(3),'on');
for m=1:Nm,plot(ax(3),Vscan,Tsol(:,m),'LineWidth',LW,'Color',C(m,:));end
xlabel(ax(3),'Common voltage (V)');
ylabel(ax(3),'{\it T}_{sol}','Interpreter','tex');
legend(ax(3),ModelName,'Location','best','Box','off','FontSize',FS1-4);
xlim(ax(3),[-10 10]);ylim(ax(3),[0 1]);

% (d) V1-only directional solar response
ax(4)=nexttile;hold(ax(4),'on');
for m=1:Nm,plot(ax(4),Vscan,DAsol(:,m),'LineWidth',LW,'Color',C(m,:));end
yline(ax(4),0,'k:','HandleVisibility','off');
xlabel(ax(4),'{\it V}_1 only (V)','Interpreter','tex');
ylabel(ax(4),'\Delta {\it A}_{sol}','Interpreter','tex');
legend(ax(4),ModelName,'Location','best','Box','off','FontSize',FS1-4);
xlim(ax(4),[-10 10]);

% (e) Vc versus EO coefficient
ax(5)=nexttile;hold(ax(5),'on');
for m=1:Nm,plot(ax(5),scaleScan,VcScale(:,m),'-o','LineWidth',LW,...
        'MarkerSize',6,'Color',C(m,:));end
yline(ax(5),0,'k:','HandleVisibility','off');
xline(ax(5),1,':','HandleVisibility','off');
xlabel(ax(5),'Electro-optic coefficient scale');
ylabel(ax(5),'Reversal voltage {\it V}_c (V)','Interpreter','tex');
legend(ax(5),ModelName,'Location','best','Box','off','FontSize',FS1-4);
xlim(ax(5),[scaleScan(1) scaleScan(end)]);ylim(ax(5),[-10 10]);

% (f) Delta Tsol sensitivity
ax(6)=nexttile;hold(ax(6),'on');
for m=1:Nm,plot(ax(6),scaleScan,DTsolScale(:,m),'-o','LineWidth',LW,...
        'MarkerSize',6,'Color',C(m,:));end
xline(ax(6),1,':','HandleVisibility','off');
xlabel(ax(6),'Electro-optic coefficient scale');
ylabel(ax(6),'\Delta {\it T}_{sol}','Interpreter','tex');
legend(ax(6),ModelName,'Location','best','Box','off','FontSize',FS1-4);
xlim(ax(6),[scaleScan(1) scaleScan(end)]);

% (g) Delta Tlum sensitivity
ax(7)=nexttile;hold(ax(7),'on');
for m=1:Nm,plot(ax(7),scaleScan,DTlumScale(:,m),'-o','LineWidth',LW,...
        'MarkerSize',6,'Color',C(m,:));end
xline(ax(7),1,':','HandleVisibility','off');
xlabel(ax(7),'Electro-optic coefficient scale');
ylabel(ax(7),'\Delta {\it T}_{lum}','Interpreter','tex');
legend(ax(7),ModelName,'Location','best','Box','off','FontSize',FS1-4);
xlim(ax(7),[scaleScan(1) scaleScan(end)]);

% (h) Directional contrast sensitivity
ax(8)=nexttile;hold(ax(8),'on');
for m=1:Nm,plot(ax(8),scaleScan,MaxDAScale(:,m),'-o','LineWidth',LW,...
        'MarkerSize',6,'Color',C(m,:));end
xline(ax(8),1,':','HandleVisibility','off');
xlabel(ax(8),'Electro-optic coefficient scale');
ylabel(ax(8),'max|\Delta {\it A}_{sol}|','Interpreter','tex');
legend(ax(8),ModelName,'Location','best','Box','off','FontSize',FS1-4);
xlim(ax(8),[scaleScan(1) scaleScan(end)]);

Panel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
for j=1:8
    set(ax(j),'FontName',FN,'FontSize',FS1,'FontWeight','normal','LineWidth',LW,...
        'TickDir','in','TickLength',[.012 .012],'Box','on','Layer','top');
    text(ax(j),0,1.06,Panel{j},'Units','normalized','FontName',FN,...
        'FontSize',FS2,'FontWeight','normal','HorizontalAlignment','left',...
        'VerticalAlignment','bottom','Clipping','off','HandleVisibility','off');
end
fprintf('Elapsed time: %.2f s\n',toc);

%% Local functions
function [Tlum,Tsol,DAsol,RecErr]=calcMetrics(wave,M,Vscan,d_nm,dDAST,ID,n0,dn_du,...
    scale,mode,anchorShift,S,Veye,iSol,iLum)
Nv=numel(Vscan);Tlum=zeros(Nv,1);Tsol=Tlum;DAsol=Tlum;RecErr=0;
for iv=1:Nv
    V=Vscan(iv);

    % Common-voltage response
    [~,Tf,~,~,Tb,~]=calcState(wave,M,[V V V],d_nm,dDAST,ID,n0,dn_du,...
        scale,mode,anchorShift);
    Tlum(iv)=weighted(Tf,wave,Veye,iLum);
    Tsol(iv)=weighted(Tf,wave,S,iSol);
    RecErr=max(RecErr,max(abs(Tf-Tb)));

    % V1-only directional response
    [~,~,Af,~,~,Ab]=calcState(wave,M,[V 0 0],d_nm,dDAST,ID,n0,dn_du,...
        scale,mode,anchorShift);
    DAsol(iv)=weighted(Af-Ab,wave,S,iSol);
end
end

function [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,V,d_nm,dDAST,ID,n0,dn_du,...
    scale,mode,anchorShift)
N=numel(wave);Rf=zeros(N,1);Tf=Rf;Af=Rf;Rb=Rf;Tb=Rf;Ab=Rf;
for iw=1:N
    d1=dastIndexModel(M(iw,ID.DAST),V(1),dDAST,n0,dn_du,scale,mode,anchorShift);
    d2=dastIndexModel(M(iw,ID.DAST),V(2),dDAST,n0,dn_du,scale,mode,anchorShift);
    d3=dastIndexModel(M(iw,ID.DAST),V(3),dDAST,n0,dn_du,scale,mode,anchorShift);

    nf=[M(iw,ID.Au),d1,M(iw,ID.ITO),M(iw,ID.SiO2),M(iw,ID.ITO),...
        d2,M(iw,ID.Ag),M(iw,ID.SiO2),M(iw,ID.Ag),d3,M(iw,ID.ITO)];

    lam=wave(iw)*1e-9;df=d_nm*1e-9;
    [Rf(iw),Tf(iw)]=tmmRT(1,1,nf,df,lam);
    [Rb(iw),Tb(iw)]=tmmRT(1,1,fliplr(nf),fliplr(df),lam);
    Af(iw)=1-Rf(iw)-Tf(iw);
    Ab(iw)=1-Rb(iw)-Tb(iw);
end
end

function nV=dastIndexModel(nk0,V,t_nm,n0,dn_du,scale,mode,anchorShift)
dn=scale*dn_du*V/t_nm;k=abs(imag(nk0));
if mode==1
    nr=n0;
elseif mode==2
    nr=real(nk0);
else
    nr=real(nk0)+anchorShift;
end
nV=(nr+dn)-1i*k;
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

function Vc=zeroCrossNearest(V,y,Vref)
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