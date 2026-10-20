% exportgraphics(gcf,'DAST_final_mechanism.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

FS1=20;FS2=26;FN='Times New Roman';LW=1.5;FigWidth=1900;FigHeight=900;mu=char(181);
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

%% Final nominal structure
% Air/Au/DAST1/ITO1/SiO2_1/ITO2/DAST2/Ag1/SiO2_2/Ag2/DAST3/ITO3/Air
dAu=5;dDAST=100;dITO=40;dSi1=90;dSi2=115;dAg=8;           % nm
d_nm=[dAu dDAST dITO dSi1 dITO dDAST dAg dSi2 dAg dDAST dITO];
df=d_nm*1e-9;
LayerLabel={'Au','DAST_1','ITO_1','SiO_2^{(1)}','ITO_2','DAST_2',...
    'Ag_1','SiO_2^{(2)}','Ag_2','DAST_3','ITO_3'};
LayerPlain={'Au','DAST1','ITO1','SiO2_1','ITO2','DAST2',...
    'Ag1','SiO2_2','Ag2','DAST3','ITO3'};
NL=numel(d_nm);

%% Sellmeier-hybrid DAST model
nConst=2.2;dn_du=3.41;
qX=1.645;nInfX=2.078;lambda0X=533;                         % nm
BlendLow=650;BlendHigh=750;
Vgrid=(-10:0.5:10)';
FieldBand=[700 2500];                                      % nm
PointsPerLayer=40;

%% Solar spectrum
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
[ws,S0]=readWeight(fileSol,3);
S=interp1(ws,S0,wave,'linear',0);
iSol=wave>=400&wave<=2500;
DenSol=trapz(wave(iSol),S(iSol));
assert(DenSol>0,'Invalid solar spectrum.');

%% Find exact Vc for V1-only control
fprintf('\nSearching final V1 reversal voltage...\n');
DAgrid=zeros(size(Vgrid));
for iv=1:numel(Vgrid)
    DAgrid(iv)=solarDirectional(Vgrid(iv),wave,M,d_nm,dDAST,ID,nConst,dn_du,...
        qX,nInfX,lambda0X,BlendLow,BlendHigh,S,iSol);
end
fun=@(v)solarDirectional(v,wave,M,d_nm,dDAST,ID,nConst,dn_du,...
    qX,nInfX,lambda0X,BlendLow,BlendHigh,S,iSol);
Vc=findExactRoot(fun,Vgrid,DAgrid);
Vstate=[-10 Vc 10];
StateText={sprintf('{\\it V}_1=-10 V'),sprintf('{\\it V}_1={\\it V}_c=%.3f V',Vc),...
    sprintf('{\\it V}_1=+10 V')};

%% Full mechanism spectra
Spec=cell(1,3);DA=zeros(numel(wave),3);Cspec=DA;Ccum=DA;
AsolF=zeros(1,3);AsolB=zeros(1,3);DAsol=zeros(1,3);
LayerSolF=zeros(3,NL);LayerSolB=zeros(3,NL);DeltaLayerSol=zeros(3,NL);
VIS=zeros(1,3);NIR=zeros(1,3);

fprintf('Calculating spectra and layer-resolved solar absorption...\n');
for s=1:3
    Spec{s}=mechanismSpectra(wave,M,[Vstate(s) 0 0],df,dDAST,ID,nConst,dn_du,...
        qX,nInfX,lambda0X,BlendLow,BlendHigh);

    DA(:,s)=Spec{s}.Af-Spec{s}.Ab;
    Cspec(:,s)=S.*DA(:,s)/DenSol;
    Ccum(:,s)=cumtrapz(wave,Cspec(:,s));

    AsolF(s)=weighted(Spec{s}.Af,wave,S,iSol);
    AsolB(s)=weighted(Spec{s}.Ab,wave,S,iSol);
    DAsol(s)=AsolF(s)-AsolB(s);

    for j=1:NL
        LayerSolF(s,j)=weighted(Spec{s}.LayerF(:,j),wave,S,iSol);
        LayerSolB(s,j)=weighted(Spec{s}.LayerB(:,j),wave,S,iSol);
    end
    DeltaLayerSol(s,:)=LayerSolF(s,:)-LayerSolB(s,:);

    ivis=wave>=400&wave<700;inir=wave>=700&wave<=2500;
    VIS(s)=trapz(wave(ivis),Cspec(ivis,s));
    NIR(s)=trapz(wave(inir),Cspec(inir,s));
end

%% Select one common mechanism wavelength
idBand=wave>=FieldBand(1)&wave<=FieldBand(2);
changeMetric=abs(Cspec(:,3)-Cspec(:,1));
tmp=changeMetric;tmp(~idBand)=-inf;
[~,iField]=max(tmp);
lambdaField=wave(iField);
fprintf('Common mechanism wavelength = %.1f nm\n',lambdaField);

%% Local fields at the same wavelength for all three states
Field=cell(1,3);
for s=1:3
    nf=makeStack(lambdaField,M(wave==lambdaField,:),[Vstate(s) 0 0],dDAST,ID,...
        nConst,dn_du,qX,nInfX,lambda0X,BlendLow,BlendHigh);
    lam=lambdaField*1e-9;

    [zF,EF]=fieldProfile(1,1,nf,df,lam,PointsPerLayer);
    [zr,Er]=fieldProfile(1,1,fliplr(nf),fliplr(df),lam,PointsPerLayer);
    zB=sum(df)-zr;[zB,ord]=sort(zB);Er=Er(ord);

    Field{s}.zF=zF*1e6;Field{s}.EF=EF;
    Field{s}.zB=zB*1e6;Field{s}.EB=Er;
end

%% Voltage-dependent weighted absorption
AsolFV=zeros(size(Vgrid));AsolBV=AsolFV;DAV=AsolFV;
for iv=1:numel(Vgrid)
    [~,~,Af,~,~,Ab]=calcState(wave,M,[Vgrid(iv) 0 0],d_nm,dDAST,ID,nConst,dn_du,...
        qX,nInfX,lambda0X,BlendLow,BlendHigh);
    AsolFV(iv)=weighted(Af,wave,S,iSol);
    AsolBV(iv)=weighted(Ab,wave,S,iSol);
    DAV(iv)=AsolFV(iv)-AsolBV(iv);
end

%% Numerical checks and layer contribution table
maxRec=0;maxClosure=0;minLayerA=inf;
for s=1:3
    maxRec=max(maxRec,max(abs(Spec{s}.Tf-Spec{s}.Tb)));
    maxClosure=max(maxClosure,max(abs(sum(Spec{s}.LayerF,2)-Spec{s}.Af)));
    maxClosure=max(maxClosure,max(abs(sum(Spec{s}.LayerB,2)-Spec{s}.Ab)));
    minLayerA=min([minLayerA;Spec{s}.LayerF(:);Spec{s}.LayerB(:)]);
end

fprintf('\n================ FINAL MECHANISM ANALYSIS =================\n');
fprintf('Structure: Au(5)/DAST(100)/ITO(40)/SiO2(90)/ITO(40)/DAST(100)/Ag(8)/SiO2(115)/Ag(8)/DAST(100)/ITO(40)\n');
fprintf('Reversal voltage Vc        = %+.4f V\n',Vc);
fprintf('Common field wavelength    = %.1f nm\n',lambdaField);
fprintf('Maximum |Tf-Tb|            = %.3e\n',maxRec);
fprintf('Layer-absorption closure   = %.3e\n',maxClosure);
fprintf('Minimum calculated Aj      = %.3e\n',minLayerA);

for s=1:3
    fprintf('\nState %d: V1=%+.4f V\n',s,Vstate(s));
    fprintf('  Asol,F / Asol,B          = %.5f / %.5f\n',AsolF(s),AsolB(s));
    fprintf('  DeltaAsol                = %+.5f\n',DAsol(s));
    fprintf('  VIS contribution         = %+.5f\n',VIS(s));
    fprintf('  NIR contribution         = %+.5f\n',NIR(s));
    fprintf('  sum_j DeltaAsol,j        = %+.5f, error=%.3e\n',...
        sum(DeltaLayerSol(s,:)),sum(DeltaLayerSol(s,:))-DAsol(s));

    [~,ord]=sort(abs(DeltaLayerSol(s,:)),'descend');
    fprintf('  Dominant directional layers:\n');
    for k=1:4
        j=ord(k);
        fprintf('    %-7s  DeltaAsol,j=%+.5f  AF=%.5f  AB=%.5f\n',...
            LayerPlain{j},DeltaLayerSol(s,j),LayerSolF(s,j),LayerSolB(s,j));
    end
end

fprintf('\nSolar-weighted layer directional contributions:\n');
fprintf('%-10s  %12s  %12s  %12s\n','Layer','V1=-10 V','V1=Vc','V1=+10 V');
for j=1:NL
    fprintf('%-10s  %+12.5f  %+12.5f  %+12.5f\n',...
        LayerPlain{j},DeltaLayerSol(1,j),DeltaLayerSol(2,j),DeltaLayerSol(3,j));
end
fprintf('===========================================================\n');

%% Figure
figure('Color','w','Position',[30 40 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
ax=gobjects(1,8);C=lines(3);zBound=[0 cumsum(d_nm)]/1000;

% (a) Spectral directional absorptance
ax(1)=nexttile;hold(ax(1),'on');
for s=1:3,plot(ax(1),wave/1000,DA(:,s),'LineWidth',LW,'Color',C(s,:));end
yline(ax(1),0,'k:','HandleVisibility','off');
xline(ax(1),lambdaField/1000,':','LineWidth',1,'HandleVisibility','off');
xlabel(ax(1),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(1),'\Delta {\it A}_{dir}={\it A}_F-{\it A}_B','Interpreter','tex');
legend(ax(1),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(1),[.4 2.5]);

% (b) Solar-weighted spectral contribution
ax(2)=nexttile;hold(ax(2),'on');
for s=1:3,plot(ax(2),wave/1000,Cspec(:,s)*1000,'LineWidth',LW,'Color',C(s,:));end
yline(ax(2),0,'k:','HandleVisibility','off');
xline(ax(2),lambdaField/1000,':','LineWidth',1,'HandleVisibility','off');
xlabel(ax(2),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(2),'Solar-weighted contribution (10^{-3} nm^{-1})');
legend(ax(2),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(2),[.4 2.5]);

% (c) Cumulative solar directional contribution
ax(3)=nexttile;hold(ax(3),'on');
for s=1:3,plot(ax(3),wave/1000,Ccum(:,s),'LineWidth',LW,'Color',C(s,:));end
yline(ax(3),0,'k:','HandleVisibility','off');
xlabel(ax(3),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(3),'Cumulative \Delta {\it A}_{sol}','Interpreter','tex');
legend(ax(3),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(3),[.4 2.5]);

% (d) Solar-weighted layer directional contribution
ax(4)=nexttile;
bar(ax(4),1:NL,DeltaLayerSol','grouped');
yline(ax(4),0,'k:','HandleVisibility','off');
set(ax(4),'XTick',1:NL,'XTickLabel',LayerLabel,'TickLabelInterpreter','tex');
xtickangle(ax(4),45);
ylabel(ax(4),'\Delta {\it A}_{{sol},j}','Interpreter','tex');
legend(ax(4),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-4);
xlim(ax(4),[.4 NL+.6]);

% (e)-(g) Same-wavelength local fields
for s=1:3
    ax(4+s)=nexttile;hold(ax(4+s),'on');
    plot(ax(4+s),Field{s}.zF,Field{s}.EF,'-','LineWidth',LW);
    plot(ax(4+s),Field{s}.zB,Field{s}.EB,'--','LineWidth',LW);
    for k=2:numel(zBound)-1
        xline(ax(4+s),zBound(k),':','Color',[.65 .65 .65],...
            'LineWidth',.7,'HandleVisibility','off');
    end
    xlabel(ax(4+s),['Position {\it z} (' mu 'm)'],'Interpreter','tex');
    ylabel(ax(4+s),'|{\it E}|^2/|{\it E}_{inc}|^2','Interpreter','tex');
    legend(ax(4+s),{'Forward','Reverse'},'Location','best','Box','off');
    xlim(ax(4+s),[0 zBound(end)]);
    text(ax(4+s),.03,.96,sprintf('%s,  \\lambda=%.3f %sm',...
        StateText{s},lambdaField/1000,mu),'Units','normalized','Interpreter','tex',...
        'FontName',FN,'FontSize',FS1-3,'VerticalAlignment','top');
end

% (h) Solar-weighted forward/reverse absorption versus V1
ax(8)=nexttile;hold(ax(8),'on');
yyaxis(ax(8),'left');
plot(ax(8),Vgrid,AsolFV,'-o','LineWidth',LW,'MarkerSize',4);
plot(ax(8),Vgrid,AsolBV,'-s','LineWidth',LW,'MarkerSize',4);
ylabel(ax(8),'{\it A}_{sol}','Interpreter','tex');

yyaxis(ax(8),'right');
plot(ax(8),Vgrid,DAV,'-^','LineWidth',LW,'MarkerSize',4);
yline(ax(8),0,'k:','HandleVisibility','off');
xline(ax(8),Vc,':','LineWidth',1.2,'HandleVisibility','off');
ylabel(ax(8),'\Delta {\it A}_{sol}','Interpreter','tex');
xlabel(ax(8),'{\it V}_1 (V)','Interpreter','tex');
legend(ax(8),{'{\it A}_{sol,F}','{\it A}_{sol,B}','\Delta {\it A}_{sol}'},...
    'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);
xlim(ax(8),[-10 10]);

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
function S=mechanismSpectra(wave,M,V,df,tDAST,ID,nConst,dn_du,q,nInf,l0,b1,b2)
N=numel(wave);NL=numel(df);
S.Rf=zeros(N,1);S.Tf=S.Rf;S.Af=S.Rf;S.Rb=S.Rf;S.Tb=S.Rf;S.Ab=S.Rf;
S.LayerF=zeros(N,NL);S.LayerB=zeros(N,NL);
for iw=1:N
    nf=makeStack(wave(iw),M(iw,:),V,tDAST,ID,nConst,dn_du,q,nInf,l0,b1,b2);
    lam=wave(iw)*1e-9;
    [S.Rf(iw),S.Tf(iw),aF]=tmmLayerAbs(1,1,nf,df,lam);
    [S.Rb(iw),S.Tb(iw),aBR]=tmmLayerAbs(1,1,fliplr(nf),fliplr(df),lam);
    S.Af(iw)=1-S.Rf(iw)-S.Tf(iw);S.Ab(iw)=1-S.Rb(iw)-S.Tb(iw);
    S.LayerF(iw,:)=aF(:).';S.LayerB(iw,:)=flipud(aBR(:)).';
end
end

function [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,V,d_nm,tDAST,ID,nConst,dn_du,q,nInf,l0,b1,b2)
N=numel(wave);Rf=zeros(N,1);Tf=Rf;Af=Rf;Rb=Rf;Tb=Rf;Ab=Rf;df=d_nm*1e-9;
for iw=1:N
    nf=makeStack(wave(iw),M(iw,:),V,tDAST,ID,nConst,dn_du,q,nInf,l0,b1,b2);
    lam=wave(iw)*1e-9;
    [Rf(iw),Tf(iw)]=tmmRT(1,1,nf,df,lam);
    [Rb(iw),Tb(iw)]=tmmRT(1,1,fliplr(nf),fliplr(df),lam);
    Af(iw)=1-Rf(iw)-Tf(iw);Ab(iw)=1-Rb(iw)-Tb(iw);
end
end

function nf=makeStack(w,row,V,tDAST,ID,nConst,dn_du,q,nInf,l0,b1,b2)
d1=dastIndex(w,row(ID.DAST),V(1),tDAST,nConst,dn_du,q,nInf,l0,b1,b2);
d2=dastIndex(w,row(ID.DAST),V(2),tDAST,nConst,dn_du,q,nInf,l0,b1,b2);
d3=dastIndex(w,row(ID.DAST),V(3),tDAST,nConst,dn_du,q,nInf,l0,b1,b2);
nf=[row(ID.Au),d1,row(ID.ITO),row(ID.SiO2),row(ID.ITO),...
    d2,row(ID.Ag),row(ID.SiO2),row(ID.Ag),d3,row(ID.ITO)];
end

function nV=dastIndex(w,nk0,V,t,nConst,dn_du,q,nInf,l0,b1,b2)
nr=dastBaseN(w,nConst,q,nInf,l0,b1,b2);
nV=(nr+dn_du*V/t)-1i*abs(imag(nk0));
end

function n=dastBaseN(w,nConst,q,nInf,l0,b1,b2)
if w<=b1,n=nConst;
elseif w>=b2,n=sellmeierX(w,q,nInf,l0);
else
    s=.5-.5*cos(pi*(w-b1)/(b2-b1));
    n=(1-s)*nConst+s*sellmeierX(w,q,nInf,l0);
end
end

function n=sellmeierX(w,q,nInf,l0)
n=sqrt(nInf^2+q*l0^2/(w^2-l0^2));
end

function val=solarDirectional(V1,wave,M,d_nm,tDAST,ID,nConst,dn_du,q,nInf,l0,b1,b2,S,iSol)
[~,~,Af,~,~,Ab]=calcState(wave,M,[V1 0 0],d_nm,tDAST,ID,nConst,dn_du,q,nInf,l0,b1,b2);
val=weighted(Af-Ab,wave,S,iSol);
end

function Vc=findExactRoot(fun,V,y)
roots=[];
for k=1:numel(V)-1
    if abs(y(k))<1e-12,roots(end+1)=V(k);end
    if y(k)*y(k+1)<0,roots(end+1)=fzero(fun,[V(k) V(k+1)]);end
end
if abs(y(end))<1e-12,roots(end+1)=V(end);end
if isempty(roots),error('No V1 directional reversal exists in [-10,+10] V.');end
[~,i]=min(abs(roots));Vc=roots(i);
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

function [R,T,a]=tmmLayerAbs(n0,ns,nf,df,lam)
N=numel(nf);MM=eye(2);Mj=cell(N,1);
for j=1:N
    n=nf(j);d=2*pi*n*df(j)/lam;
    Mj{j}=[cos(d) 1i*sin(d)/n;1i*n*sin(d) cos(d)];
    MM=MM*Mj{j};
end
den=n0*MM(1,1)+n0*ns*MM(1,2)+MM(2,1)+ns*MM(2,2);
r=(n0*MM(1,1)+n0*ns*MM(1,2)-MM(2,1)-ns*MM(2,2))/den;
t=2*n0/den;R=real(abs(r)^2);T=real(max(real(ns),0)/max(real(n0),eps)*abs(t)^2);
v=[1+r;n0*(1-r)];a=zeros(N,1);
for j=1:N
    Pin=pflux(v,n0);v2=Mj{j}\v;Pout=pflux(v2,n0);
    a(j)=real(Pin-Pout);v=v2;
end
a(abs(a)<1e-12)=0;
end

function [z,E2]=fieldProfile(n0,ns,nf,df,lam,Np)
N=numel(nf);MM=eye(2);
for j=1:N
    n=nf(j);d=2*pi*n*df(j)/lam;
    MM=MM*[cos(d) 1i*sin(d)/n;1i*n*sin(d) cos(d)];
end
den=n0*MM(1,1)+n0*ns*MM(1,2)+MM(2,1)+ns*MM(2,2);
r=(n0*MM(1,1)+n0*ns*MM(1,2)-MM(2,1)-ns*MM(2,2))/den;
v=[1+r;n0*(1-r)];z=[];E2=[];z0=0;
for j=1:N
    n=nf(j);x=linspace(0,df(j),Np+1);if j>1,x=x(2:end);end
    for k=1:numel(x)
        d=2*pi*n*x(k)/lam;
        Mx=[cos(d) 1i*sin(d)/n;1i*n*sin(d) cos(d)];
        vx=Mx\v;
        z(end+1,1)=z0+x(k);E2(end+1,1)=abs(vx(1))^2;
    end
    d=2*pi*n*df(j)/lam;
    v=[cos(d) 1i*sin(d)/n;1i*n*sin(d) cos(d)]\v;
    z0=z0+df(j);
end
end

function p=pflux(v,nInc)
p=real(v(1)*conj(v(2)))/max(real(nInc),eps);
end

function y=weighted(x,w,W,idx)
x=x(:);w=w(:);W=W(:);
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