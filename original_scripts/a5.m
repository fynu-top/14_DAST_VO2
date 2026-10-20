% exportgraphics(gcf,'DAST_mechanism.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

%% Figure controls
FS1=20;FS2=26;FN='Times New Roman';LW=1.5;
FigWidth=1900;FigHeight=900;mu=char(181);
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

%% Fixed final structure
% Air/Au/DAST1/ITO1/SiO2_1/ITO2/DAST2/Ag1/SiO2_2/Ag2/DAST3/ITO3/Air
dAu=5;dDAST=100;dITO=40;dSi1=80;dSi2=90;dAg=8;             % nm
d_nm=[dAu dDAST dITO dSi1 dITO dDAST dAg dSi2 dAg dDAST dITO];
df=d_nm*1e-9;
n0DAST=2.2;dn_du=3.41;                                     % nm/V

LayerLabel={'Au','DAST_1','ITO_1','SiO_2^{(1)}','ITO_2','DAST_2',...
            'Ag_1','SiO_2^{(2)}','Ag_2','DAST_3','ITO_3'};
LayerPlain={'Au','DAST1','ITO1','SiO2_1','ITO2','DAST2',...
            'Ag1','SiO2_2','Ag2','DAST3','ITO3'};
NL=numel(LayerLabel);

%% Solar weighting
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
[ws,S0]=readWeight(fileSol,3);
S=interp1(ws,S0,wave,'linear',0);
iSol=wave>=400&wave<=2500;
assert(trapz(wave(iSol),S(iSol))>0,'Invalid solar weighting data.');

%% Automatically determine V1 reversal voltage
Vcoarse=(-10:1:10)';
DAcoarse=zeros(size(Vcoarse));
fprintf('\nSearching for V1 directional-absorption reversal...\n');
for k=1:numel(Vcoarse)
    DAcoarse(k)=solarDirAtV(Vcoarse(k),wave,M,d_nm,dDAST,ID,n0DAST,dn_du,S,iSol);
end

ib=find(DAcoarse(1:end-1).*DAcoarse(2:end)<0,1);
if isempty(ib)
    iz=find(abs(DAcoarse)<1e-12,1);
    if isempty(iz),error('No V1 reversal found over -10 to +10 V.');else,Vc=Vcoarse(iz);end
else
    fun=@(v)solarDirAtV(v,wave,M,d_nm,dDAST,ID,n0DAST,dn_du,S,iSol);
    Vc=fzero(fun,[Vcoarse(ib) Vcoarse(ib+1)]);
end

fprintf('Calculated V1 reversal voltage: Vc = %.4f V\n',Vc);

%% Smooth DeltaAsol(V1) curve
Vplot=(-10:0.5:10)';
DAplot=zeros(size(Vplot));
for k=1:numel(Vplot)
    DAplot(k)=solarDirAtV(Vplot(k),wave,M,d_nm,dDAST,ID,n0DAST,dn_du,S,iSol);
end

%% Three mechanism states
Vstate=[-10 Vc 10];
StateText={sprintf('{\\it V}_1=-10 V'),...
           sprintf('{\\it V}_1={\\it V}_c=%.2f V',Vc),...
           sprintf('{\\it V}_1=+10 V')};

Spec=cell(1,3);
LambdaRep=zeros(1,3);
DAsol=zeros(1,3);
AFLayer=zeros(3,NL);
ABLayer=zeros(3,NL);
Field=cell(1,3);

fprintf('\nCalculating layer absorption and local fields...\n');
for s=1:3
    V1=Vstate(s);
    Spec{s}=mechanismSpectra(wave,M,[V1 0 0],df,dDAST,ID,n0DAST,dn_du);

    DA=Spec{s}.Af-Spec{s}.Ab;
    [~,iRep]=max(abs(DA));
    LambdaRep(s)=wave(iRep);

    DAsol(s)=weighted(DA,wave,S,iSol);
    AFLayer(s,:)=Spec{s}.LayerF(iRep,:);
    ABLayer(s,:)=Spec{s}.LayerB(iRep,:);

    nf=makeStack(M(iRep,:),[V1 0 0],dDAST,ID,n0DAST,dn_du);
    lam=wave(iRep)*1e-9;

    [zF,EF]=fieldProfile(1,1,nf,df,lam,35);
    [zr,Er]=fieldProfile(1,1,fliplr(nf),fliplr(df),lam,35);

    zB=sum(df)-zr;
    [zB,ord]=sort(zB);Er=Er(ord);

    Field{s}.zF=zF/1e-6;
    Field{s}.EF=EF;
    Field{s}.zB=zB/1e-6;
    Field{s}.EB=Er;

    sumF=sum(AFLayer(s,:));
    sumB=sum(ABLayer(s,:));
    AtF=Spec{s}.Af(iRep);
    AtB=Spec{s}.Ab(iRep);

    [~,iLF]=max(AFLayer(s,:));
    [~,iLB]=max(ABLayer(s,:));

    fprintf('\nState %d: V1 = %+.4f V\n',s,V1);
    fprintf('  lambda*                   = %.1f nm\n',LambdaRep(s));
    fprintf('  DeltaAdir(lambda*)        = %+.6f\n',DA(iRep));
    fprintf('  DeltaAsol                 = %+.6f\n',DAsol(s));
    fprintf('  Forward: R/T/A            = %.6f / %.6f / %.6f\n',...
        Spec{s}.Rf(iRep),Spec{s}.Tf(iRep),AtF);
    fprintf('  Reverse: R/T/A            = %.6f / %.6f / %.6f\n',...
        Spec{s}.Rb(iRep),Spec{s}.Tb(iRep),AtB);
    fprintf('  Sum layer A, forward      = %.6f, error = %.3e\n',sumF,sumF-AtF);
    fprintf('  Sum layer A, reverse      = %.6f, error = %.3e\n',sumB,sumB-AtB);
    fprintf('  Dominant forward layer    = %s (%.5f)\n',LayerPlain{iLF},AFLayer(s,iLF));
    fprintf('  Dominant reverse layer    = %s (%.5f)\n',LayerPlain{iLB},ABLayer(s,iLB));
end

%% Global numerical checks
maxRec=0;maxClosure=0;minLayerA=inf;
for s=1:3
    maxRec=max(maxRec,max(abs(Spec{s}.Tf-Spec{s}.Tb)));
    maxClosure=max(maxClosure,max(abs(sum(Spec{s}.LayerF,2)-Spec{s}.Af)));
    maxClosure=max(maxClosure,max(abs(sum(Spec{s}.LayerB,2)-Spec{s}.Ab)));
    minLayerA=min([minLayerA;Spec{s}.LayerF(:);Spec{s}.LayerB(:)]);
end

fprintf('\n================ MECHANISM CHECK =================\n');
fprintf('Vc                       = %.4f V\n',Vc);
fprintf('Representative lambda    = %.1f / %.1f / %.1f nm\n',LambdaRep);
fprintf('DeltaAsol                = %+.5f / %+.5f / %+.5f\n',DAsol);
fprintf('max |Tf-Tb|              = %.3e\n',maxRec);
fprintf('max layer closure error  = %.3e\n',maxClosure);
fprintf('minimum layer A          = %.3e\n',minLayerA);
fprintf('==================================================\n');

%% Figure
figure('Color','w','Position',[30 40 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
ax=gobjects(1,8);C=lines(3);CF=[0 .447 .741];CB=[.850 .325 .098];

% (a) Spectral directional absorptance
ax(1)=nexttile;hold(ax(1),'on');
for s=1:3
    DA=Spec{s}.Af-Spec{s}.Ab;
    plot(ax(1),wave/1000,DA,'LineWidth',LW,'Color',C(s,:));
    [~,ii]=min(abs(wave-LambdaRep(s)));
    plot(ax(1),LambdaRep(s)/1000,DA(ii),'o','Color',C(s,:),...
        'MarkerFaceColor','w','MarkerSize',7,'LineWidth',LW,'HandleVisibility','off');
end
yline(ax(1),0,'k:','HandleVisibility','off');
xlabel(ax(1),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(1),'\Delta {\it A}_{dir}={\it A}_F-{\it A}_B','Interpreter','tex');
legend(ax(1),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(1),[.4 2.5]);

% (b) Solar-weighted directional reversal
ax(2)=nexttile;hold(ax(2),'on');
plot(ax(2),Vplot,DAplot,'-','LineWidth',LW);
yline(ax(2),0,'k:','HandleVisibility','off');
xline(ax(2),Vc,':','LineWidth',LW,'HandleVisibility','off');
for s=1:3
    plot(ax(2),Vstate(s),DAsol(s),'o','MarkerFaceColor','w','MarkerSize',7,'LineWidth',LW);
end
xlabel(ax(2),'{\it V}_1 (V)','Interpreter','tex');
ylabel(ax(2),'\Delta {\it A}_{sol}','Interpreter','tex');
xlim(ax(2),[-10 10]);

% (c)-(e) Layer-resolved absorption
for s=1:3
    ax(2+s)=nexttile;hold(ax(2+s),'on');
    plot(ax(2+s),1:NL,AFLayer(s,:),'-o','Color',CF,'LineWidth',LW,'MarkerSize',5);
    plot(ax(2+s),1:NL,ABLayer(s,:),'--s','Color',CB,'LineWidth',LW,'MarkerSize',5);
    set(ax(2+s),'XTick',1:NL,'XTickLabel',LayerLabel,'TickLabelInterpreter','tex');
    xtickangle(ax(2+s),45);
    ylabel(ax(2+s),'Layer absorptance {\it A}_j','Interpreter','tex');
    xlim(ax(2+s),[.5 NL+.5]);
    legend(ax(2+s),{'Forward','Reverse'},'Location','best','Box','off','FontSize',FS1-4);
    text(ax(2+s),.03,.96,sprintf('%s,  \\lambda^*=%.3f %sm',...
        StateText{s},LambdaRep(s)/1000,mu),'Units','normalized','Interpreter','tex',...
        'FontName',FN,'FontSize',FS1-3,'VerticalAlignment','top');
end

% (f)-(h) Local field
zBound=[0 cumsum(d_nm)]/1000;                          % um
for s=1:3
    ax(5+s)=nexttile;hold(ax(5+s),'on');
    plot(ax(5+s),Field{s}.zF,Field{s}.EF,'-','Color',CF,'LineWidth',LW);
    plot(ax(5+s),Field{s}.zB,Field{s}.EB,'--','Color',CB,'LineWidth',LW);

    for k=2:numel(zBound)-1
        xline(ax(5+s),zBound(k),':','Color',[.55 .55 .55],...
            'LineWidth',.8,'HandleVisibility','off');
    end

    xlabel(ax(5+s),['Position {\it z} (' mu 'm)'],'Interpreter','tex');
    ylabel(ax(5+s),'|{\it E}|^2/|{\it E}_{inc}|^2','Interpreter','tex');
    xlim(ax(5+s),[0 zBound(end)]);
    legend(ax(5+s),{'Forward','Reverse'},'Location','best','Box','off','FontSize',FS1-4);

    text(ax(5+s),.03,.96,sprintf('%s,  \\lambda^*=%.3f %sm',...
        StateText{s},LambdaRep(s)/1000,mu),'Units','normalized','Interpreter','tex',...
        'FontName',FN,'FontSize',FS1-3,'VerticalAlignment','top');
end

Panel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
for j=1:8
    set(ax(j),'FontName',FN,'FontSize',FS1,'FontWeight','normal',...
        'LineWidth',LW,'TickDir','in','TickLength',[.012 .012],...
        'Box','on','Layer','top');
    text(ax(j),0,1.06,Panel{j},'Units','normalized','FontName',FN,...
        'FontSize',FS2,'FontWeight','normal','HorizontalAlignment','left',...
        'VerticalAlignment','bottom','Clipping','off','HandleVisibility','off');
end

fprintf('Elapsed time: %.2f s\n',toc);

%% Local functions
function S=mechanismSpectra(wave,M,V,df,tDAST,ID,n0,dn_du)
N=numel(wave);NL=numel(df);
S.Rf=zeros(N,1);S.Tf=S.Rf;S.Af=S.Rf;S.Rb=S.Rf;S.Tb=S.Rf;S.Ab=S.Rf;
S.LayerF=zeros(N,NL);S.LayerB=zeros(N,NL);

for iw=1:N
    nf=makeStack(M(iw,:),V,tDAST,ID,n0,dn_du);
    lam=wave(iw)*1e-9;

    [S.Rf(iw),S.Tf(iw),aF]=tmmLayerAbs(1,1,nf,df,lam);
    [S.Rb(iw),S.Tb(iw),aBR]=tmmLayerAbs(1,1,fliplr(nf),fliplr(df),lam);

    S.Af(iw)=1-S.Rf(iw)-S.Tf(iw);
    S.Ab(iw)=1-S.Rb(iw)-S.Tb(iw);

    S.LayerF(iw,:)=aF(:).';
    S.LayerB(iw,:)=flipud(aBR(:)).';
end
end

function nf=makeStack(row,V,tDAST,ID,n0,dn_du)
d1=dastIndex(row(ID.DAST),V(1),tDAST,n0,dn_du);
d2=dastIndex(row(ID.DAST),V(2),tDAST,n0,dn_du);
d3=dastIndex(row(ID.DAST),V(3),tDAST,n0,dn_du);

nf=[row(ID.Au),d1,row(ID.ITO),row(ID.SiO2),row(ID.ITO),...
    d2,row(ID.Ag),row(ID.SiO2),row(ID.Ag),d3,row(ID.ITO)];
end

function nV=dastIndex(nk0,V,t_nm,n0,dn_du)
nV=(n0+dn_du*V/t_nm)-1i*abs(imag(nk0));
end

function val=solarDirAtV(V1,wave,M,d_nm,tDAST,ID,n0,dn_du,S,iSol)
N=numel(wave);DA=zeros(N,1);df=d_nm*1e-9;
for iw=1:N
    nf=makeStack(M(iw,:),[V1 0 0],tDAST,ID,n0,dn_du);
    lam=wave(iw)*1e-9;

    [Rf,Tf]=tmmRT(1,1,nf,df,lam);
    [Rb,Tb]=tmmRT(1,1,fliplr(nf),fliplr(df),lam);

    DA(iw)=(1-Rf-Tf)-(1-Rb-Tb);
end
val=weighted(DA,wave,S,iSol);
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

function [R,T,a]=tmmLayerAbs(n0,ns,nf,df,lam)
N=numel(nf);M=eye(2);Mj=cell(N,1);

for j=1:N
    n=nf(j);delta=2*pi*n*df(j)/lam;
    Mj{j}=[cos(delta) 1i*sin(delta)/n;1i*n*sin(delta) cos(delta)];
    M=M*Mj{j};
end

den=n0*M(1,1)+n0*ns*M(1,2)+M(2,1)+ns*M(2,2);
r=(n0*M(1,1)+n0*ns*M(1,2)-M(2,1)-ns*M(2,2))/den;
t=2*n0/den;

R=real(abs(r)^2);
T=real(max(real(ns),0)/max(real(n0),eps)*abs(t)^2);

v=[1+r;n0*(1-r)];
a=zeros(N,1);

for j=1:N
    Pin=pflux(v,n0);
    v2=Mj{j}\v;
    Pout=pflux(v2,n0);
    a(j)=real(Pin-Pout);
    v=v2;
end

a(abs(a)<1e-11)=0;
end

function [z,E2]=fieldProfile(n0,ns,nf,df,lam,Np)
N=numel(nf);M=eye(2);

for j=1:N
    n=nf(j);delta=2*pi*n*df(j)/lam;
    M=M*[cos(delta) 1i*sin(delta)/n;1i*n*sin(delta) cos(delta)];
end

den=n0*M(1,1)+n0*ns*M(1,2)+M(2,1)+ns*M(2,2);
r=(n0*M(1,1)+n0*ns*M(1,2)-M(2,1)-ns*M(2,2))/den;

v=[1+r;n0*(1-r)];
z=[];E2=[];z0=0;

for j=1:N
    n=nf(j);
    x=linspace(0,df(j),Np+1);

    if j>1,x=x(2:end);end

    for k=1:numel(x)
        delta=2*pi*n*x(k)/lam;
        Mx=[cos(delta) 1i*sin(delta)/n;1i*n*sin(delta) cos(delta)];
        vx=Mx\v;

        z(end+1,1)=z0+x(k);
        E2(end+1,1)=abs(vx(1))^2;
    end

    delta=2*pi*n*df(j)/lam;
    Mj=[cos(delta) 1i*sin(delta)/n;1i*n*sin(delta) cos(delta)];
    v=Mj\v;
    z0=z0+df(j);
end
end

function p=pflux(v,n0)
p=real(v(1)*conj(v(2)))/max(real(n0),eps);
end

function y=weighted(x,w,W,idx)
x=x(:);w=w(:);W=W(:);
y=trapz(w(idx),x(idx).*W(idx))/trapz(w(idx),W(idx));
end

function file=findData(root,names)
project=fileparts(root);
folders={root,project,fullfile(root,'Atmospheric data'),...
    fullfile(project,'Atmospheric data')};
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