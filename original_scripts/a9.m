% exportgraphics(gcf,'SiO2_Sellmeier_reversal_optimization.pdf','ContentType','vector');
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
assert(size(M,2)>=13,'nk_data must contain Ag as material ID 13.');

%% Fixed structure
% Air/Au/DAST1/ITO/SiO2_1/ITO/DAST2/Ag/SiO2_2/Ag/DAST3/ITO/Air
dAu=5;dAg=8;dITO=40;dDAST=100;                             % nm
baseSi1=80;baseSi2=90;

%% Sellmeier-hybrid DAST
nConst=2.2;dn_du=3.41;
qX=1.645;nInfX=2.078;lambda0X=533;                         % nm
BlendLow=650;BlendHigh=750;

%% Scan controls
Si1Scan=60:5:120;
Si2Scan=80:5:140;
Vscan=(-10:1:10)';
N1=numel(Si1Scan);N2=numel(Si2Scan);

% Smart-window/reversal constraints
MinTlum0=0.55;
MinDTsol=0.018;
MinDAminus=0.015;
MaxDAplus=-0.005;
MaxAbsVc=7;

%% Spectral weights
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
assert(trapz(wave(iLum),Veye(iLum))>0&&trapz(wave(iSol),S(iSol))>0,'Invalid spectral weights.');

%% Allocate maps
VcMap=nan(N2,N1);Tlum0Map=zeros(N2,N1);DTsolMap=zeros(N2,N1);
DAminusMap=zeros(N2,N1);DAplusMap=zeros(N2,N1);MarginMap=zeros(N2,N1);
MaxDAMap=zeros(N2,N1);RecMap=zeros(N2,N1);

fprintf('\nScanning Sellmeier-hybrid SiO2 spacer designs...\n');
for i1=1:N1
    for i2=1:N2
        s1=Si1Scan(i1);s2=Si2Scan(i2);
        d_nm=[dAu dDAST dITO s1 dITO dDAST dAg s2 dAg dDAST dITO];

        Tlum=zeros(numel(Vscan),1);Tsol=Tlum;DA=Tlum;rec=0;
        for iv=1:numel(Vscan)
            V=Vscan(iv);

            % Common-voltage smart-window response
            [~,Tf,~,~,Tb,~]=calcState(wave,M,[V V V],d_nm,dDAST,ID,nConst,dn_du,...
                qX,nInfX,lambda0X,BlendLow,BlendHigh);
            Tlum(iv)=weighted(Tf,wave,Veye,iLum);
            Tsol(iv)=weighted(Tf,wave,S,iSol);
            rec=max(rec,max(abs(Tf-Tb)));

            % V1-only directional solar response
            [~,~,Af,~,~,Ab]=calcState(wave,M,[V 0 0],d_nm,dDAST,ID,nConst,dn_du,...
                qX,nInfX,lambda0X,BlendLow,BlendHigh);
            DA(iv)=weighted(Af-Ab,wave,S,iSol);
        end

        VcMap(i2,i1)=zeroCrossNearest(Vscan,DA,0);
        Tlum0Map(i2,i1)=Tlum(Vscan==0);
        DTsolMap(i2,i1)=max(Tsol)-min(Tsol);
        DAminusMap(i2,i1)=DA(Vscan==-10);
        DAplusMap(i2,i1)=DA(Vscan==10);
        MarginMap(i2,i1)=min(DAminusMap(i2,i1),-DAplusMap(i2,i1));
        MaxDAMap(i2,i1)=max(abs(DA));
        RecMap(i2,i1)=rec;
    end
    fprintf('SiO2_1 = %3g nm completed\n',Si1Scan(i1));
end

%% Feasible designs
Feasible=Tlum0Map>=MinTlum0 & DTsolMap>=MinDTsol & ...
         DAminusMap>=MinDAminus & DAplusMap<=MaxDAplus & ...
         abs(VcMap)<=MaxAbsVc & isfinite(VcMap);

fprintf('\n================ SELLMЕIER SPACER OPTIMIZATION =================\n');
fprintf('Total designs        : %d\n',N1*N2);
fprintf('Feasible designs     : %d\n',sum(Feasible,'all'));
fprintf('Maximum reciprocity error = %.3e\n',max(RecMap,[],'all'));

%% Baseline
[~,ib1]=min(abs(Si1Scan-baseSi1));[~,ib2]=min(abs(Si2Scan-baseSi2));

%% Select three useful candidates
valid=find(Feasible);

if isempty(valid)
    warning('No design satisfies every hard constraint. Selecting from reversal-capable designs.');
    Candidate=isfinite(VcMap)&MarginMap>0&Tlum0Map>=0.50;
    valid=find(Candidate);
end

if isempty(valid),error('No directional-reversal design found in the scanned range.');end

% 1) strongest reversal
[~,k]=max(MarginMap(valid));iStrong=valid(k);

% 2) lowest-voltage reversal
[~,k]=min(abs(VcMap(valid)));iLowV=valid(k);

% 3) balanced point: equal normalized preference for margin, low |Vc|, Tlum0, DeltaTsol
m=MarginMap(valid);v=abs(VcMap(valid));tl=Tlum0Map(valid);dt=DTsolMap(valid);
mn=normalize01(m);vn=1-normalize01(v);tn=normalize01(tl);dn=normalize01(dt);
Score=0.35*mn+0.30*vn+0.20*tn+0.15*dn;
[~,k]=max(Score);iBalanced=valid(k);

[rStrong,cStrong]=ind2sub(size(VcMap),iStrong);
[rLow,cLow]=ind2sub(size(VcMap),iLowV);
[rBal,cBal]=ind2sub(size(VcMap),iBalanced);

fprintf('\nBaseline:\n');
printDesign(baseSi1,baseSi2,VcMap(ib2,ib1),Tlum0Map(ib2,ib1),DTsolMap(ib2,ib1),...
    DAminusMap(ib2,ib1),DAplusMap(ib2,ib1),MarginMap(ib2,ib1));

fprintf('\nStrong-reversal candidate:\n');
printDesign(Si1Scan(cStrong),Si2Scan(rStrong),VcMap(rStrong,cStrong),...
    Tlum0Map(rStrong,cStrong),DTsolMap(rStrong,cStrong),...
    DAminusMap(rStrong,cStrong),DAplusMap(rStrong,cStrong),MarginMap(rStrong,cStrong));

fprintf('\nLow-voltage candidate:\n');
printDesign(Si1Scan(cLow),Si2Scan(rLow),VcMap(rLow,cLow),...
    Tlum0Map(rLow,cLow),DTsolMap(rLow,cLow),...
    DAminusMap(rLow,cLow),DAplusMap(rLow,cLow),MarginMap(rLow,cLow));

fprintf('\nBalanced candidate:\n');
printDesign(Si1Scan(cBal),Si2Scan(rBal),VcMap(rBal,cBal),...
    Tlum0Map(rBal,cBal),DTsolMap(rBal,cBal),...
    DAminusMap(rBal,cBal),DAplusMap(rBal,cBal),MarginMap(rBal,cBal));
fprintf('==============================================================\n');

%% Full validation: baseline vs balanced candidate
DesignSi=[baseSi1 baseSi2;Si1Scan(cBal) Si2Scan(rBal)];
Name={'Baseline','Balanced'};
Nd=2;Nv=numel(Vscan);
TlumFull=zeros(Nv,Nd);TsolFull=TlumFull;DAFull=TlumFull;
VcFull=nan(1,Nd);

for d=1:Nd
    d_nm=[dAu dDAST dITO DesignSi(d,1) dITO dDAST dAg DesignSi(d,2) dAg dDAST dITO];

    for iv=1:Nv
        V=Vscan(iv);

        [~,Tf,~,~,~,~]=calcState(wave,M,[V V V],d_nm,dDAST,ID,nConst,dn_du,...
            qX,nInfX,lambda0X,BlendLow,BlendHigh);
        TlumFull(iv,d)=weighted(Tf,wave,Veye,iLum);
        TsolFull(iv,d)=weighted(Tf,wave,S,iSol);

        [~,~,Af,~,~,Ab]=calcState(wave,M,[V 0 0],d_nm,dDAST,ID,nConst,dn_du,...
            qX,nInfX,lambda0X,BlendLow,BlendHigh);
        DAFull(iv,d)=weighted(Af-Ab,wave,S,iSol);
    end

    VcFull(d)=zeroCrossNearest(Vscan,DAFull(:,d),0);
end

%% Figure
figure('Color','w','Position',[30 50 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');ax=gobjects(1,8);

% (a) Reversal voltage
ax(1)=nexttile;
imagesc(ax(1),Si1Scan,Si2Scan,VcMap);axis(ax(1),'xy');caxis(ax(1),[-10 10]);
xlabel(ax(1),'SiO_2 spacer 1 (nm)');ylabel(ax(1),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(1));cb.Label.String='Reversal voltage {\it V}_c (V)';cb.Label.Interpreter='tex';
hold(ax(1),'on');markDesigns(ax(1),baseSi1,baseSi2,Si1Scan(cBal),Si2Scan(rBal));

% (b) Reversal margin
ax(2)=nexttile;
imagesc(ax(2),Si1Scan,Si2Scan,MarginMap);axis(ax(2),'xy');
xlabel(ax(2),'SiO_2 spacer 1 (nm)');ylabel(ax(2),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(2));cb.Label.String='Reversal margin {\it M}_{rev}';cb.Label.Interpreter='tex';
hold(ax(2),'on');markDesigns(ax(2),baseSi1,baseSi2,Si1Scan(cBal),Si2Scan(rBal));

% (c) DeltaAsol at -10 V
ax(3)=nexttile;
imagesc(ax(3),Si1Scan,Si2Scan,DAminusMap);axis(ax(3),'xy');
xlabel(ax(3),'SiO_2 spacer 1 (nm)');ylabel(ax(3),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(3));cb.Label.String='\Delta {\it A}_{sol}(-10 V)';cb.Label.Interpreter='tex';
hold(ax(3),'on');markDesigns(ax(3),baseSi1,baseSi2,Si1Scan(cBal),Si2Scan(rBal));

% (d) DeltaAsol at +10 V
ax(4)=nexttile;
imagesc(ax(4),Si1Scan,Si2Scan,DAplusMap);axis(ax(4),'xy');
xlabel(ax(4),'SiO_2 spacer 1 (nm)');ylabel(ax(4),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(4));cb.Label.String='\Delta {\it A}_{sol}(+10 V)';cb.Label.Interpreter='tex';
hold(ax(4),'on');markDesigns(ax(4),baseSi1,baseSi2,Si1Scan(cBal),Si2Scan(rBal));

% (e) Tlum at zero bias
ax(5)=nexttile;
imagesc(ax(5),Si1Scan,Si2Scan,Tlum0Map);axis(ax(5),'xy');caxis(ax(5),[0 1]);
xlabel(ax(5),'SiO_2 spacer 1 (nm)');ylabel(ax(5),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(5));cb.Label.String='{\it T}_{lum}(0 V)';cb.Label.Interpreter='tex';
hold(ax(5),'on');markDesigns(ax(5),baseSi1,baseSi2,Si1Scan(cBal),Si2Scan(rBal));

% (f) DeltaTsol
ax(6)=nexttile;
imagesc(ax(6),Si1Scan,Si2Scan,DTsolMap);axis(ax(6),'xy');
xlabel(ax(6),'SiO_2 spacer 1 (nm)');ylabel(ax(6),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(6));cb.Label.String='\Delta {\it T}_{sol}';cb.Label.Interpreter='tex';
hold(ax(6),'on');markDesigns(ax(6),baseSi1,baseSi2,Si1Scan(cBal),Si2Scan(rBal));

% (g) Trade-off plane
ax(7)=nexttile;hold(ax(7),'on');
ok=isfinite(VcMap)&MarginMap>0;
scatter(ax(7),abs(VcMap(ok)),MarginMap(ok),38,Tlum0Map(ok),'filled');
plot(ax(7),abs(VcMap(ib2,ib1)),MarginMap(ib2,ib1),'ks','MarkerSize',9,'LineWidth',1.6);
plot(ax(7),abs(VcMap(rBal,cBal)),MarginMap(rBal,cBal),'kp','MarkerSize',12,'LineWidth',1.8);
xlabel(ax(7),'|{\it V}_c| (V)','Interpreter','tex');
ylabel(ax(7),'Reversal margin {\it M}_{rev}','Interpreter','tex');
cb=colorbar(ax(7));cb.Label.String='{\it T}_{lum}(0 V)';cb.Label.Interpreter='tex';
legend(ax(7),{'Reversal designs','Baseline','Balanced'},'Location','best','Box','off');

% (h) Full validation
ax(8)=nexttile;hold(ax(8),'on');
yyaxis(ax(8),'left');
plot(ax(8),Vscan,DAFull(:,1),'--','LineWidth',LW);
plot(ax(8),Vscan,DAFull(:,2),'-','LineWidth',LW);
yline(ax(8),0,'k:','HandleVisibility','off');
ylabel(ax(8),'\Delta {\it A}_{sol}','Interpreter','tex');

yyaxis(ax(8),'right');
plot(ax(8),Vscan,TsolFull(:,1),'--','LineWidth',LW);
plot(ax(8),Vscan,TsolFull(:,2),'-','LineWidth',LW);
ylabel(ax(8),'{\it T}_{sol}','Interpreter','tex');
xlabel(ax(8),'Voltage (V)');
legend(ax(8),{'Baseline \Delta{\it A}_{sol}','Balanced \Delta{\it A}_{sol}',...
    'Baseline {\it T}_{sol}','Balanced {\it T}_{sol}'},...
    'Location','best','Box','off','Interpreter','tex','FontSize',FS1-4);
xlim(ax(8),[-10 10]);

Panel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
for j=1:8
    set(ax(j),'FontName',FN,'FontSize',FS1,'FontWeight','normal','LineWidth',LW,...
        'TickDir','in','TickLength',[.012 .012],'Box','on','Layer','top');
    text(ax(j),0,1.06,Panel{j},'Units','normalized','FontName',FN,'FontSize',FS2,...
        'FontWeight','normal','HorizontalAlignment','left','VerticalAlignment','bottom',...
        'Clipping','off','HandleVisibility','off');
end

fprintf('\nFull validation:\n');
for d=1:Nd
    fprintf('%s: Si1=%g nm, Si2=%g nm, Vc=%+.3f V, Tlum0=%.4f, DeltaTsol=%.4f, DA(-10)=%+.4f, DA(+10)=%+.4f\n',...
        Name{d},DesignSi(d,1),DesignSi(d,2),VcFull(d),TlumFull(Vscan==0,d),...
        max(TsolFull(:,d))-min(TsolFull(:,d)),DAFull(1,d),DAFull(end,d));
end
fprintf('Elapsed time: %.2f s\n',toc);

%% Local functions
function [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,V,d_nm,tDAST,ID,nConst,dn_du,q,nInf,l0,b1,b2)
N=numel(wave);Rf=zeros(N,1);Tf=Rf;Af=Rf;Rb=Rf;Tb=Rf;Ab=Rf;
for iw=1:N
    d1=dastIndex(wave(iw),M(iw,ID.DAST),V(1),tDAST,nConst,dn_du,q,nInf,l0,b1,b2);
    d2=dastIndex(wave(iw),M(iw,ID.DAST),V(2),tDAST,nConst,dn_du,q,nInf,l0,b1,b2);
    d3=dastIndex(wave(iw),M(iw,ID.DAST),V(3),tDAST,nConst,dn_du,q,nInf,l0,b1,b2);

    nf=[M(iw,ID.Au),d1,M(iw,ID.ITO),M(iw,ID.SiO2),M(iw,ID.ITO),...
        d2,M(iw,ID.Ag),M(iw,ID.SiO2),M(iw,ID.Ag),d3,M(iw,ID.ITO)];

    lam=wave(iw)*1e-9;df=d_nm*1e-9;
    [Rf(iw),Tf(iw)]=tmmRT(1,1,nf,df,lam);
    [Rb(iw),Tb(iw)]=tmmRT(1,1,fliplr(nf),fliplr(df),lam);
    Af(iw)=1-Rf(iw)-Tf(iw);Ab(iw)=1-Rb(iw)-Tb(iw);
end
end

function nV=dastIndex(w,nk0,V,t,nConst,dn_du,q,nInf,l0,b1,b2)
nr=dastBaseN(w,nConst,q,nInf,l0,b1,b2);
nV=(nr+dn_du*V/t)-1i*abs(imag(nk0));
end

function n=dastBaseN(w,nConst,q,nInf,l0,b1,b2)
if w<=b1
    n=nConst;
elseif w>=b2
    n=sellmeierX(w,q,nInf,l0);
else
    s=.5-.5*cos(pi*(w-b1)/(b2-b1));
    n=(1-s)*nConst+s*sellmeierX(w,q,nInf,l0);
end
end

function n=sellmeierX(w,q,nInf,l0)
n=sqrt(nInf^2+q*l0^2/(w^2-l0^2));
end

function [R,T]=tmmRT(n0,ns,nf,df,lam)
M=eye(2);
for j=1:numel(nf)
    n=nf(j);d=2*pi*n*df(j)/lam;
    M=M*[cos(d) 1i*sin(d)/n;1i*n*sin(d) cos(d)];
end
den=n0*M(1,1)+n0*ns*M(1,2)+M(2,1)+ns*M(2,2);
r=(n0*M(1,1)+n0*ns*M(1,2)-M(2,1)-ns*M(2,2))/den;
t=2*n0/den;
R=real(abs(r)^2);T=real(max(real(ns),0)/max(real(n0),eps)*abs(t)^2);
end

function y=weighted(x,w,W,idx)
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

function y=normalize01(x)
if max(x)-min(x)<eps,y=ones(size(x));else,y=(x-min(x))/(max(x)-min(x));end
end

function printDesign(s1,s2,Vc,Tlum0,DTsol,DAm,DAp,Margin)
fprintf('  SiO2_1 / SiO2_2 = %g / %g nm\n',s1,s2);
fprintf('  Vc                = %+.3f V\n',Vc);
fprintf('  Tlum(0 V)         = %.4f\n',Tlum0);
fprintf('  DeltaTsol         = %.4f\n',DTsol);
fprintf('  DeltaAsol(-10 V)  = %+.4f\n',DAm);
fprintf('  DeltaAsol(+10 V)  = %+.4f\n',DAp);
fprintf('  Reversal margin   = %.4f\n',Margin);
end

function markDesigns(ax,b1,b2,s1,s2)
plot(ax,b1,b2,'ks','MarkerSize',9,'LineWidth',1.5);
plot(ax,s1,s2,'kp','MarkerSize',12,'LineWidth',1.8);
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