% exportgraphics(gcf,'DAST_TE_TM_angle.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

%% Global controls
FS1=20;FS2=26;FN='Times New Roman';LW=1.5;FigWidth=1900;FigHeight=900;
thetaScan=(0:5:75)';Vscan=(-10:1:10)';xGap=zeros(1,8);yGap=zeros(1,8);
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

%% Sellmeier-hybrid DAST model
nConst=2.2;dn_du=3.41;qX=1.645;nInfX=2.078;lambda0X=533;
BlendLow=650;BlendHigh=750;
nBase=zeros(size(wave));
for iw=1:numel(wave)
    nBase(iw)=dastBaseN(wave(iw),nConst,qX,nInfX,lambda0X,BlendLow,BlendHigh);
end
kDAST=abs(imag(M(:,ID.DAST)));

%% Spectral weights
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
assert(trapz(wave(iLum),Veye(iLum))>0&&trapz(wave(iSol),S(iSol))>0,'Invalid spectral weights.');

%% Exact normal-incidence Vc
fprintf('\nSearching exact normal-incidence reversal voltage...\n');
DA0=zeros(size(Vscan));
for iv=1:numel(Vscan)
    DA0(iv)=solarDirectionalAngle(Vscan(iv),0,'TE',wave,M,d_nm,dDAST,ID,nBase,kDAST,dn_du,S,iSol);
end
fun=@(v)solarDirectionalAngle(v,0,'TE',wave,M,d_nm,dDAST,ID,nBase,kDAST,dn_du,S,iSol);
Vc0=findExactRoot(fun,Vscan,DA0,0);
Vstate=[-10 Vc0 10];Ns=3;
StateText={sprintf('{\\it V}_1=-10 V'),sprintf('{\\it V}_1={\\it V}_c=%.3f V',Vc0),...
    sprintf('{\\it V}_1=+10 V')};
fprintf('Normal-incidence Vc = %+.4f V\n',Vc0);

%% Angular response
Pol={'TE','TM'};Np=2;Na=numel(thetaScan);
Tlum=zeros(Na,Ns,Np);Tsol=Tlum;AsF=Tlum;AsB=Tlum;DA=Tlum;
RecState=zeros(Na,Ns,Np);VcAngle=nan(Na,Np);RecScan=zeros(Na,Np);

fprintf('\nCalculating TE/TM angular response...\n');
for p=1:Np
    fprintf('%s polarization...\n',Pol{p});
    Vref=Vc0;
    for ia=1:Na
        th=thetaScan(ia);

        % Three fixed voltage states
        for s=1:Ns
            [~,Tf,Af,~,Tb,Ab]=calcStateAngle(wave,M,[Vstate(s) 0 0],d_nm,dDAST,ID,...
                nBase,kDAST,dn_du,th,Pol{p});
            Tlum(ia,s,p)=weighted(Tf,wave,Veye,iLum);
            Tsol(ia,s,p)=weighted(Tf,wave,S,iSol);
            AsF(ia,s,p)=weighted(Af,wave,S,iSol);
            AsB(ia,s,p)=weighted(Ab,wave,S,iSol);
            DA(ia,s,p)=AsF(ia,s,p)-AsB(ia,s,p);
            RecState(ia,s,p)=max(abs(Tf-Tb));
        end

        % Angle-dependent V1 reversal voltage
        daV=zeros(size(Vscan));rr=0;
        for iv=1:numel(Vscan)
            [~,Tf,Af,~,Tb,Ab]=calcStateAngle(wave,M,[Vscan(iv) 0 0],d_nm,dDAST,ID,...
                nBase,kDAST,dn_du,th,Pol{p});
            daV(iv)=weighted(Af-Ab,wave,S,iSol);
            rr=max(rr,max(abs(Tf-Tb)));
        end
        roots=rootCandidates(Vscan,daV);
        if ~isempty(roots)
            [~,ir]=min(abs(roots-Vref));VcAngle(ia,p)=roots(ir);Vref=VcAngle(ia,p);
        end
        RecScan(ia,p)=rr;
    end
end

%% Angular robustness metrics
RevMask=false(Na,Np);
for p=1:Np
    RevMask(:,p)=DA(:,1,p)>0 & DA(:,3,p)<0;
end

MaxStableAngle=zeros(1,Np);
for p=1:Np
    bad=find(~RevMask(:,p),1);
    if isempty(bad),MaxStableAngle(p)=thetaScan(end);
    elseif bad==1,MaxStableAngle(p)=NaN;
    else,MaxStableAngle(p)=thetaScan(bad-1);
    end
end

maxRec=max([RecState(:);RecScan(:)]);
maxPolTlum=max(abs(Tlum(:,:,1)-Tlum(:,:,2)),[],'all');
maxPolTsol=max(abs(Tsol(:,:,1)-Tsol(:,:,2)),[],'all');
maxPolDA=max(abs(DA(:,:,1)-DA(:,:,2)),[],'all');

fprintf('\n================ ANGULAR PERFORMANCE =================\n');
fprintf('Scalar DAST model: Sellmeier-hybrid n_x(lambda,V)\n');
fprintf('Normal-incidence Vc             = %+.4f V\n',Vc0);
fprintf('Maximum reciprocity error       = %.3e\n',maxRec);
fprintf('Maximum TE/TM difference Tlum   = %.4f\n',maxPolTlum);
fprintf('Maximum TE/TM difference Tsol   = %.4f\n',maxPolTsol);
fprintf('Maximum TE/TM difference DAsol  = %.4f\n',maxPolDA);
fprintf('Reversal retained continuously:\n');
fprintf('  TE: 0-%g deg\n',MaxStableAngle(1));
fprintf('  TM: 0-%g deg\n',MaxStableAngle(2));

checkAng=[0 30 45 60 75];
fprintf('\n Angle     Vc_TE      Vc_TM      DA_TE(-10)  DA_TE(+10)  DA_TM(-10)  DA_TM(+10)\n');
for q=1:numel(checkAng)
    ia=find(thetaScan==checkAng(q),1);
    fprintf('%5.0f   %9s   %9s     %+8.4f    %+8.4f    %+8.4f    %+8.4f\n',...
        thetaScan(ia),numText(VcAngle(ia,1)),numText(VcAngle(ia,2)),...
        DA(ia,1,1),DA(ia,3,1),DA(ia,1,2),DA(ia,3,2));
end
fprintf('======================================================\n');

%% Figure
figure('Color','w','Position',[30 40 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
ax=gobjects(1,8);C=lines(3);

% (a) TE luminous transmittance
ax(1)=nexttile;hold(ax(1),'on');
for s=1:Ns,plot(ax(1),thetaScan,Tlum(:,s,1),'LineWidth',LW,'Color',C(s,:));end
xlabel(ax(1),'Incidence angle (deg)');ylabel(ax(1),'{\it T}_{lum}^{TE}','Interpreter','tex');
lg=legend(ax(1),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);shiftLegend(lg,xGap(1),yGap(1));
xlim(ax(1),[0 75]);ylim(ax(1),[0 1]);

% (b) TM luminous transmittance
ax(2)=nexttile;hold(ax(2),'on');
for s=1:Ns,plot(ax(2),thetaScan,Tlum(:,s,2),'LineWidth',LW,'Color',C(s,:));end
xlabel(ax(2),'Incidence angle (deg)');ylabel(ax(2),'{\it T}_{lum}^{TM}','Interpreter','tex');
lg=legend(ax(2),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);shiftLegend(lg,xGap(2),yGap(2));
xlim(ax(2),[0 75]);ylim(ax(2),[0 1]);

% (c) TE solar transmittance
ax(3)=nexttile;hold(ax(3),'on');
for s=1:Ns,plot(ax(3),thetaScan,Tsol(:,s,1),'LineWidth',LW,'Color',C(s,:));end
xlabel(ax(3),'Incidence angle (deg)');ylabel(ax(3),'{\it T}_{sol}^{TE}','Interpreter','tex');
lg=legend(ax(3),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);shiftLegend(lg,xGap(3),yGap(3));
xlim(ax(3),[0 75]);ylim(ax(3),[0 1]);

% (d) TM solar transmittance
ax(4)=nexttile;hold(ax(4),'on');
for s=1:Ns,plot(ax(4),thetaScan,Tsol(:,s,2),'LineWidth',LW,'Color',C(s,:));end
xlabel(ax(4),'Incidence angle (deg)');ylabel(ax(4),'{\it T}_{sol}^{TM}','Interpreter','tex');
lg=legend(ax(4),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);shiftLegend(lg,xGap(4),yGap(4));
xlim(ax(4),[0 75]);ylim(ax(4),[0 1]);

% (e) TE forward/reverse solar absorptance
ax(5)=nexttile;hold(ax(5),'on');lab={};
for s=1:Ns
    plot(ax(5),thetaScan,AsF(:,s,1),'-','LineWidth',LW,'Color',C(s,:));
    plot(ax(5),thetaScan,AsB(:,s,1),'--','LineWidth',LW,'Color',C(s,:));
    lab{end+1}=sprintf('%s, F',StateText{s});lab{end+1}=sprintf('%s, B',StateText{s});
end
xlabel(ax(5),'Incidence angle (deg)');ylabel(ax(5),'{\it A}_{sol}^{TE}','Interpreter','tex');
lg=legend(ax(5),lab,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-5);shiftLegend(lg,xGap(5),yGap(5));
xlim(ax(5),[0 75]);ylim(ax(5),[0 1]);

% (f) TM forward/reverse solar absorptance
ax(6)=nexttile;hold(ax(6),'on');lab={};
for s=1:Ns
    plot(ax(6),thetaScan,AsF(:,s,2),'-','LineWidth',LW,'Color',C(s,:));
    plot(ax(6),thetaScan,AsB(:,s,2),'--','LineWidth',LW,'Color',C(s,:));
    lab{end+1}=sprintf('%s, F',StateText{s});lab{end+1}=sprintf('%s, B',StateText{s});
end
xlabel(ax(6),'Incidence angle (deg)');ylabel(ax(6),'{\it A}_{sol}^{TM}','Interpreter','tex');
lg=legend(ax(6),lab,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-5);shiftLegend(lg,xGap(6),yGap(6));
xlim(ax(6),[0 75]);ylim(ax(6),[0 1]);

% (g) Directional solar absorptance difference: TE/TM
ax(7)=nexttile;hold(ax(7),'on');lab={};
for s=1:Ns
    plot(ax(7),thetaScan,DA(:,s,1),'-','LineWidth',LW,'Color',C(s,:));
    plot(ax(7),thetaScan,DA(:,s,2),'--','LineWidth',LW,'Color',C(s,:));
    lab{end+1}=sprintf('%s, TE',StateText{s});lab{end+1}=sprintf('%s, TM',StateText{s});
end
yline(ax(7),0,'k:','HandleVisibility','off');
xlabel(ax(7),'Incidence angle (deg)');ylabel(ax(7),'\Delta {\it A}_{sol}','Interpreter','tex');
lg=legend(ax(7),lab,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-5);shiftLegend(lg,xGap(7),yGap(7));
xlim(ax(7),[0 75]);

% (h) Angle-dependent reversal voltage
ax(8)=nexttile;hold(ax(8),'on');
plot(ax(8),thetaScan,VcAngle(:,1),'-o','LineWidth',LW,'MarkerSize',5);
plot(ax(8),thetaScan,VcAngle(:,2),'-s','LineWidth',LW,'MarkerSize',5);
yline(ax(8),0,'k:','HandleVisibility','off');
yline(ax(8),Vc0,'k--','LineWidth',1,'HandleVisibility','off');
xlabel(ax(8),'Incidence angle (deg)');ylabel(ax(8),'Reversal voltage {\it V}_c (V)','Interpreter','tex');
lg=legend(ax(8),{'TE','TM'},'Location','best','Box','off');shiftLegend(lg,xGap(8),yGap(8));
xlim(ax(8),[0 75]);ylim(ax(8),[-10 10]);

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
function [Rf,Tf,Af,Rb,Tb,Ab]=calcStateAngle(wave,M,V,d_nm,tD,ID,nBase,kDAST,dn_du,theta,pol)
N=numel(wave);Rf=zeros(N,1);Tf=Rf;Af=Rf;Rb=Rf;Tb=Rf;Ab=Rf;df=d_nm*1e-9;
for iw=1:N
    d1=(nBase(iw)+dn_du*V(1)/tD)-1i*kDAST(iw);
    d2=(nBase(iw)+dn_du*V(2)/tD)-1i*kDAST(iw);
    d3=(nBase(iw)+dn_du*V(3)/tD)-1i*kDAST(iw);
    nf=[M(iw,ID.Au),d1,M(iw,ID.ITO),M(iw,ID.SiO2),M(iw,ID.ITO),...
        d2,M(iw,ID.Ag),M(iw,ID.SiO2),M(iw,ID.Ag),d3,M(iw,ID.ITO)];
    lam=wave(iw)*1e-9;
    [Rf(iw),Tf(iw)]=tmmRTAngle(1,1,nf,df,lam,theta,pol);
    [Rb(iw),Tb(iw)]=tmmRTAngle(1,1,fliplr(nf),fliplr(df),lam,theta,pol);
    Af(iw)=real(1-Rf(iw)-Tf(iw));Ab(iw)=real(1-Rb(iw)-Tb(iw));
end
end

function val=solarDirectionalAngle(V1,theta,pol,wave,M,d_nm,tD,ID,nBase,kDAST,dn_du,S,iSol)
[~,~,Af,~,~,Ab]=calcStateAngle(wave,M,[V1 0 0],d_nm,tD,ID,nBase,kDAST,dn_du,theta,pol);
val=weighted(Af-Ab,wave,S,iSol);
end

function [R,T]=tmmRTAngle(n0,ns,nf,df,lam,theta0,pol)
s0=n0*sind(theta0);c0=cosd(theta0);cs=layerCos(ns,s0);
eta0=admittance(n0,c0,pol);etas=admittance(ns,cs,pol);
MM=eye(2);
for j=1:numel(nf)
    n=nf(j);c=layerCos(n,s0);eta=admittance(n,c,pol);
    delta=2*pi*n*c*df(j)/lam;
    MM=MM*[cos(delta) 1i*sin(delta)/eta;1i*eta*sin(delta) cos(delta)];
end
B=MM(1,1)+MM(1,2)*etas;C=MM(2,1)+MM(2,2)*etas;
den=eta0*B+C;r=(eta0*B-C)/den;t=2*eta0/den;
R=real(abs(r)^2);
T=real(real(etas)/real(eta0)*abs(t)^2);
end

function c=layerCos(n,s0)
c=sqrt(1-(s0./n).^2);
if real(n*c)<0,c=-c;end
end

function eta=admittance(n,c,pol)
if strcmpi(pol,'TE'),eta=n*c;
else,eta=n/c;
end
end

function n=dastBaseN(w,nConst,q,nInf,l0,b1,b2)
if w<=b1,n=nConst;
elseif w>=b2,n=sqrt(nInf^2+q*l0^2/(w^2-l0^2));
else
    s=.5-.5*cos(pi*(w-b1)/(b2-b1));
    ns=sqrt(nInf^2+q*l0^2/(w^2-l0^2));
    n=(1-s)*nConst+s*ns;
end
end

function roots=rootCandidates(V,y)
roots=[];
for k=1:numel(V)-1
    if abs(y(k))<1e-12,roots(end+1)=V(k);end
    if y(k)*y(k+1)<0
        roots(end+1)=V(k)-y(k)*(V(k+1)-V(k))/(y(k+1)-y(k));
    end
end
if abs(y(end))<1e-12,roots(end+1)=V(end);end
roots=unique(round(roots,10));
end

function Vc=findExactRoot(fun,V,y,Vref)
roots=[];
for k=1:numel(V)-1
    if abs(y(k))<1e-12,roots(end+1)=V(k);end
    if y(k)*y(k+1)<0,roots(end+1)=fzero(fun,[V(k) V(k+1)]);end
end
if abs(y(end))<1e-12,roots(end+1)=V(end);end
if isempty(roots),error('No normal-incidence reversal found in [-10,+10] V.');end
[~,i]=min(abs(roots-Vref));Vc=roots(i);
end

function y=weighted(x,w,W,idx)
x=x(:);w=w(:);W=W(:);
y=trapz(w(idx),x(idx).*W(idx))/trapz(w(idx),W(idx));
end

function s=numText(x)
if isfinite(x),s=sprintf('%+.3f',x);else,s='   --   ';end
end

function shiftLegend(lg,dx,dy)
if isempty(lg)||~isvalid(lg),return;end
drawnow;lg.Units='normalized';p=lg.Position;p(1)=p(1)+dx;p(2)=p(2)+dy;lg.Position=p;
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
A=readmatrix(file);assert(size(A,2)>=col,'Insufficient columns in %s.',file);
ok=isfinite(A(:,1))&isfinite(A(:,col))&A(:,1)>0;
w=A(ok,1);y=A(ok,col);if max(w)<100,w=w*1000;end
[w,id]=sort(w);y=y(id);[w,id]=unique(w,'stable');y=y(id);
end