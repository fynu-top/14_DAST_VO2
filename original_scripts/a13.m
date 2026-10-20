% exportgraphics(gcf,'DAST_Vtheta_unpolarized.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

%% Global controls
FS1=20;FS2=26;FN='Times New Roman';LW=1.5;FigWidth=1900;FigHeight=900;
thetaScan=(0:2.5:75)';Vscan=(-10:0.5:10)';StrongMargin=0.005;
xGap1=0;yGap1=0;xGap2=0;yGap2=0;xGap3=0;yGap3=0;xGap4=0;yGap4=0;
xGap5=0;yGap5=0;xGap6=0;yGap6=0;xGap7=0;yGap7=0;xGap8=0;yGap8=0;
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

%% Final fixed structure
% Air/Au/DAST1/ITO1/SiO2_1/ITO2/DAST2/Ag1/SiO2_2/Ag2/DAST3/ITO3/Air
dAu=5;dDAST=100;dITO=40;dSi1=90;dSi2=115;dAg=8;           % nm
d_nm=[dAu dDAST dITO dSi1 dITO dDAST dAg dSi2 dAg dDAST dITO];

%% Sellmeier-hybrid DAST
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

%% Exact normal-incidence reversal voltage
Vroot=(-10:1:10)';DAroot=zeros(size(Vroot));
fprintf('\nSearching exact normal-incidence reversal voltage...\n');
for iv=1:numel(Vroot)
    DAroot(iv)=solarDirectionalAngle(Vroot(iv),0,'TE',wave,M,d_nm,dDAST,ID,nBase,kDAST,dn_du,S,iSol);
end
fun=@(v)solarDirectionalAngle(v,0,'TE',wave,M,d_nm,dDAST,ID,nBase,kDAST,dn_du,S,iSol);
Vc0=findExactRoot(fun,Vroot,DAroot,0);
Vstate=[-10 Vc0 10];Ns=3;
StateText={sprintf('{\\it V}_1=-10 V'),sprintf('{\\it V}_1={\\it V}_c=%.3f V',Vc0),...
    sprintf('{\\it V}_1=+10 V')};
fprintf('Normal-incidence Vc = %+.4f V\n',Vc0);

%% Full V-theta phase maps
Na=numel(thetaScan);Nv=numel(Vscan);
TlumTE=zeros(Na,Nv);TsolTE=TlumTE;AsFTE=TlumTE;AsBTE=TlumTE;DATE=TlumTE;
TlumTM=TlumTE;TsolTM=TlumTE;AsFTM=TlumTE;AsBTM=TlumTE;DATM=TlumTE;
Rec=zeros(Na,Nv,2);

fprintf('\nCalculating full V-theta phase maps...\n');
for ia=1:Na
    th=thetaScan(ia);
    for iv=1:Nv
        V=Vscan(iv);

        [~,Tf,Af,~,Tb,Ab]=calcStateAngle(wave,M,[V 0 0],d_nm,dDAST,ID,nBase,kDAST,dn_du,th,'TE');
        TlumTE(ia,iv)=weighted(Tf,wave,Veye,iLum);TsolTE(ia,iv)=weighted(Tf,wave,S,iSol);
        AsFTE(ia,iv)=weighted(Af,wave,S,iSol);AsBTE(ia,iv)=weighted(Ab,wave,S,iSol);
        DATE(ia,iv)=AsFTE(ia,iv)-AsBTE(ia,iv);Rec(ia,iv,1)=max(abs(Tf-Tb));

        [~,Tf,Af,~,Tb,Ab]=calcStateAngle(wave,M,[V 0 0],d_nm,dDAST,ID,nBase,kDAST,dn_du,th,'TM');
        TlumTM(ia,iv)=weighted(Tf,wave,Veye,iLum);TsolTM(ia,iv)=weighted(Tf,wave,S,iSol);
        AsFTM(ia,iv)=weighted(Af,wave,S,iSol);AsBTM(ia,iv)=weighted(Ab,wave,S,iSol);
        DATM(ia,iv)=AsFTM(ia,iv)-AsBTM(ia,iv);Rec(ia,iv,2)=max(abs(Tf-Tb));
    end
    fprintf('Angle = %4.1f deg completed\n',th);
end

%% Unpolarized maps
TlumU=(TlumTE+TlumTM)/2;TsolU=(TsolTE+TsolTM)/2;
AsFU=(AsFTE+AsFTM)/2;AsBU=(AsBTE+AsBTM)/2;DAU=(DATE+DATM)/2;

iVm=find(abs(Vscan+10)<1e-12,1);iVp=find(abs(Vscan-10)<1e-12,1);
MarginTE=min(DATE(:,iVm),-DATE(:,iVp));
MarginTM=min(DATM(:,iVm),-DATM(:,iVp));
MarginU=min(DAU(:,iVm),-DAU(:,iVp));

RevTE=DATE(:,iVm)>0&DATE(:,iVp)<0;
RevTM=DATM(:,iVm)>0&DATM(:,iVp)<0;
RevU =DAU(:,iVm)>0&DAU(:,iVp)<0;

StrongTE=MarginTE>=StrongMargin;StrongTM=MarginTM>=StrongMargin;StrongU=MarginU>=StrongMargin;
angRevTE=continuousAngle(thetaScan,RevTE);angRevTM=continuousAngle(thetaScan,RevTM);angRevU=continuousAngle(thetaScan,RevU);
angStrongTE=continuousAngle(thetaScan,StrongTE);angStrongTM=continuousAngle(thetaScan,StrongTM);angStrongU=continuousAngle(thetaScan,StrongU);

%% Exact three-state unpolarized angular response
TlumState=zeros(Na,Ns);TsolState=TlumState;DAState=TlumState;
for ia=1:Na
    for s=1:Ns
        [~,TfE,AfE,~,~,AbE]=calcStateAngle(wave,M,[Vstate(s) 0 0],d_nm,dDAST,ID,nBase,kDAST,dn_du,thetaScan(ia),'TE');
        [~,TfM,AfM,~,~,AbM]=calcStateAngle(wave,M,[Vstate(s) 0 0],d_nm,dDAST,ID,nBase,kDAST,dn_du,thetaScan(ia),'TM');
        TlumState(ia,s)=.5*(weighted(TfE,wave,Veye,iLum)+weighted(TfM,wave,Veye,iLum));
        TsolState(ia,s)=.5*(weighted(TfE,wave,S,iSol)+weighted(TfM,wave,S,iSol));
        DAState(ia,s)=.5*(weighted(AfE-AbE,wave,S,iSol)+weighted(AfM-AbM,wave,S,iSol));
    end
end

%% Extract all zero-crossing branches
RootsTE=cell(Na,1);RootsTM=cell(Na,1);RootsU=cell(Na,1);
for ia=1:Na
    RootsTE{ia}=allRoots(Vscan,DATE(ia,:));
    RootsTM{ia}=allRoots(Vscan,DATM(ia,:));
    RootsU{ia}=allRoots(Vscan,DAU(ia,:));
end

%% Summary
maxRec=max(Rec,[],'all');
fprintf('\n================ V-THETA / UNPOLARIZED SUMMARY ================\n');
fprintf('Normal-incidence Vc                   = %+.4f V\n',Vc0);
fprintf('Maximum reciprocity error             = %.3e\n',maxRec);
fprintf('\nPolarity reversal retained continuously:\n');
fprintf('  TE           : 0 - %.1f deg\n',angRevTE);
fprintf('  TM           : 0 - %.1f deg\n',angRevTM);
fprintf('  Unpolarized  : 0 - %.1f deg\n',angRevU);
fprintf('\nStrong reversal Mrev >= %.3f retained continuously:\n',StrongMargin);
fprintf('  TE           : 0 - %.1f deg\n',angStrongTE);
fprintf('  TM           : 0 - %.1f deg\n',angStrongTM);
fprintf('  Unpolarized  : 0 - %.1f deg\n',angStrongU);

fprintf('\n Angle   DAu(-10)   DAu(Vc0)   DAu(+10)   Mrev,u     zero roots, unpolarized\n');
checkAng=[0 15 30 35 40 45 50 55 60 65 75];
for q=1:numel(checkAng)
    [~,ia]=min(abs(thetaScan-checkAng(q)));
    fprintf('%5.1f   %+8.4f   %+8.4f   %+8.4f   %+8.4f    %s\n',...
        thetaScan(ia),DAState(ia,1),DAState(ia,2),DAState(ia,3),MarginU(ia),rootText(RootsU{ia}));
end
fprintf('================================================================\n');

%% Figure
figure('Color','w','Position',[25 40 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
ax=gobjects(1,8);
cmax=max(abs([DATE(:);DATM(:);DAU(:)]));C=lines(3);

% (a) TE V-theta phase map
ax(1)=nexttile;
imagesc(ax(1),thetaScan,Vscan,DATE');axis(ax(1),'xy');caxis(ax(1),[-cmax cmax]);hold(ax(1),'on');
contour(ax(1),thetaScan,Vscan,DATE',[0 0],'k','LineWidth',1.6);
xlabel(ax(1),'Incidence angle (deg)');ylabel(ax(1),'{\it V}_1 (V)','Interpreter','tex');
cb=colorbar(ax(1));cb.Label.String='\Delta {\it A}_{sol}^{TE}';cb.Label.Interpreter='tex';

% (b) TM V-theta phase map
ax(2)=nexttile;
imagesc(ax(2),thetaScan,Vscan,DATM');axis(ax(2),'xy');caxis(ax(2),[-cmax cmax]);hold(ax(2),'on');
contour(ax(2),thetaScan,Vscan,DATM',[0 0],'k','LineWidth',1.6);
xlabel(ax(2),'Incidence angle (deg)');ylabel(ax(2),'{\it V}_1 (V)','Interpreter','tex');
cb=colorbar(ax(2));cb.Label.String='\Delta {\it A}_{sol}^{TM}';cb.Label.Interpreter='tex';

% (c) Unpolarized V-theta phase map
ax(3)=nexttile;
imagesc(ax(3),thetaScan,Vscan,DAU');axis(ax(3),'xy');caxis(ax(3),[-cmax cmax]);hold(ax(3),'on');
contour(ax(3),thetaScan,Vscan,DAU',[0 0],'k','LineWidth',1.8);
xlabel(ax(3),'Incidence angle (deg)');ylabel(ax(3),'{\it V}_1 (V)','Interpreter','tex');
cb=colorbar(ax(3));cb.Label.String='\Delta {\it A}_{sol}^{unpol}';cb.Label.Interpreter='tex';

% (d) Reversal margin
ax(4)=nexttile;hold(ax(4),'on');
plot(ax(4),thetaScan,MarginTE,'-','LineWidth',LW);
plot(ax(4),thetaScan,MarginTM,'--','LineWidth',LW);
plot(ax(4),thetaScan,MarginU,'-o','LineWidth',LW,'MarkerSize',4);
yline(ax(4),0,'k:','HandleVisibility','off');
yline(ax(4),StrongMargin,'k--','LineWidth',1,'HandleVisibility','off');
xlabel(ax(4),'Incidence angle (deg)');
ylabel(ax(4),'Reversal margin {\it M}_{rev}','Interpreter','tex');
lg=legend(ax(4),{'TE','TM','Unpolarized'},'Location','best','Box','off');shiftLegend(lg,xGap4,yGap4);
xlim(ax(4),[0 75]);

% (e) Unpolarized luminous transmittance
ax(5)=nexttile;hold(ax(5),'on');
for s=1:Ns,plot(ax(5),thetaScan,TlumState(:,s),'LineWidth',LW,'Color',C(s,:));end
xlabel(ax(5),'Incidence angle (deg)');
ylabel(ax(5),'{\it T}_{lum}^{unpol}','Interpreter','tex');
lg=legend(ax(5),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);shiftLegend(lg,xGap5,yGap5);
xlim(ax(5),[0 75]);ylim(ax(5),[0 1]);

% (f) Unpolarized solar transmittance
ax(6)=nexttile;hold(ax(6),'on');
for s=1:Ns,plot(ax(6),thetaScan,TsolState(:,s),'LineWidth',LW,'Color',C(s,:));end
xlabel(ax(6),'Incidence angle (deg)');
ylabel(ax(6),'{\it T}_{sol}^{unpol}','Interpreter','tex');
lg=legend(ax(6),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);shiftLegend(lg,xGap6,yGap6);
xlim(ax(6),[0 75]);ylim(ax(6),[0 1]);

% (g) Unpolarized directional absorptance
ax(7)=nexttile;hold(ax(7),'on');
for s=1:Ns,plot(ax(7),thetaScan,DAState(:,s),'LineWidth',LW,'Color',C(s,:));end
yline(ax(7),0,'k:','HandleVisibility','off');
xlabel(ax(7),'Incidence angle (deg)');
ylabel(ax(7),'\Delta {\it A}_{sol}^{unpol}','Interpreter','tex');
lg=legend(ax(7),StateText,'Location','best','Box','off','Interpreter','tex','FontSize',FS1-3);shiftLegend(lg,xGap7,yGap7);
xlim(ax(7),[0 75]);

% (h) All zero-crossing branches
ax(8)=nexttile;hold(ax(8),'on');
firstTE=true;firstTM=true;firstU=true;
for ia=1:Na
    if ~isempty(RootsTE{ia})
        h=scatter(ax(8),thetaScan(ia)*ones(size(RootsTE{ia})),RootsTE{ia},28,'o');
        if firstTE,h.DisplayName='TE';firstTE=false;else,h.HandleVisibility='off';end
    end
    if ~isempty(RootsTM{ia})
        h=scatter(ax(8),thetaScan(ia)*ones(size(RootsTM{ia})),RootsTM{ia},32,'s');
        if firstTM,h.DisplayName='TM';firstTM=false;else,h.HandleVisibility='off';end
    end
    if ~isempty(RootsU{ia})
        h=scatter(ax(8),thetaScan(ia)*ones(size(RootsU{ia})),RootsU{ia},36,'^','filled');
        if firstU,h.DisplayName='Unpolarized';firstU=false;else,h.HandleVisibility='off';end
    end
end
yline(ax(8),0,'k:','HandleVisibility','off');
xlabel(ax(8),'Incidence angle (deg)');
ylabel(ax(8),'All reversal voltages {\it V}_c (V)','Interpreter','tex');
lg=legend(ax(8),'Location','best','Box','off');shiftLegend(lg,xGap8,yGap8);
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
s0=n0*sind(theta0);kz0=physicalKz(n0,s0);kzs=physicalKz(ns,s0);
c0=kz0/n0;cs=kzs/ns;eta0=admittance(n0,c0,pol);etas=admittance(ns,cs,pol);
MM=eye(2);
for j=1:numel(nf)
    n=nf(j);kz=physicalKz(n,s0);c=kz/n;eta=admittance(n,c,pol);
    delta=2*pi*kz*df(j)/lam;
    MM=MM*[cos(delta) 1i*sin(delta)/eta;1i*eta*sin(delta) cos(delta)];
end
B=MM(1,1)+MM(1,2)*etas;C=MM(2,1)+MM(2,2)*etas;
den=eta0*B+C;r=(eta0*B-C)/den;t=2*eta0/den;
R=real(abs(r)^2);T=real(real(etas)/real(eta0)*abs(t)^2);
end

function kz=physicalKz(n,s0)
kz=sqrt(n.^2-s0.^2);
if real(kz)<0||(abs(real(kz))<1e-12&&imag(kz)>0),kz=-kz;end
end

function eta=admittance(n,c,pol)
if strcmpi(pol,'TE'),eta=n*c;else,eta=n/c;end
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

function r=allRoots(V,y)
r=[];
for k=1:numel(V)-1
    if abs(y(k))<1e-10,r(end+1)=V(k);end
    if y(k)*y(k+1)<0
        r(end+1)=V(k)-y(k)*(V(k+1)-V(k))/(y(k+1)-y(k));
    end
end
if abs(y(end))<1e-10,r(end+1)=V(end);end
if ~isempty(r),r=unique(round(r,8));end
end

function Vc=findExactRoot(fun,V,y,Vref)
r=[];br=[];
for k=1:numel(V)-1
    if abs(y(k))<1e-12,r(end+1)=V(k);br(end+1,:)=V(k)*[1 1];end
    if y(k)*y(k+1)<0
        rr=V(k)-y(k)*(V(k+1)-V(k))/(y(k+1)-y(k));
        r(end+1)=rr;br(end+1,:)=[V(k) V(k+1)];
    end
end
if isempty(r),error('No reversal found in [-10,+10] V.');end
[~,i]=min(abs(r-Vref));
if br(i,1)==br(i,2),Vc=br(i,1);else,Vc=fzero(fun,br(i,:));end
end

function a=continuousAngle(theta,mask)
bad=find(~mask,1);
if isempty(bad),a=theta(end);
elseif bad==1,a=NaN;
else,a=theta(bad-1);
end
end

function s=rootText(r)
if isempty(r),s='--';return;end
s=sprintf('%+.2f, ',r);s=s(1:end-2);
end

function y=weighted(x,w,W,idx)
x=x(:);w=w(:);W=W(:);
y=trapz(w(idx),x(idx).*W(idx))/trapz(w(idx),W(idx));
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