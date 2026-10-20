% exportgraphics(gcf,'DAST_final_robustness.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

%% Global controls
FS1=20;FS2=26;FN='Times New Roman';LW=1.5;FigWidth=1900;FigHeight=880;
Nmc=300;Vstep=2;rng(1);                                    % set Nmc=500,Vstep=1 for final SI if desired
Vscan=(-10:Vstep:10)';
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

%% Final nominal design
dAu0=5;dAg0=8;dITO=40;dSi10=90;dSi20=115;dD0=100;dn_du=3.41;

%% Sellmeier-hybrid DAST model
nConst=2.2;qX=1.645;nInfX=2.078;lambda0X=533;BlendLow=650;BlendHigh=750;
nBase=zeros(size(wave));
for iw=1:numel(wave),nBase(iw)=dastBaseN(wave(iw),nConst,qX,nInfX,lambda0X,BlendLow,BlendHigh);end
kDAST=abs(imag(M(:,ID.DAST)));

%% Spectral weights
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;

%% Nominal reference
P0=[dAu0 dAg0 dAg0 dSi10 dSi20 dD0 dD0 dD0 1];
Nom=evaluateDesign(P0,wave,M,Vscan,ID,nBase,kDAST,dITO,dn_du,S,Veye,iSol,iLum);

% Refine nominal Vc with fzero
fun=@(v)solarDirectional(v,P0,wave,M,ID,nBase,kDAST,dITO,dn_du,S,iSol);
VcExact=refineVc(fun,Vscan,Nom.DA);
fprintf('\n================ NOMINAL STRUCTURE =================\n');
fprintf('Au=%.1f, DAST=%.1f, ITO=%.1f, SiO2_1=%.1f, SiO2_2=%.1f, Ag=%.1f nm\n',...
    dAu0,dD0,dITO,dSi10,dSi20,dAg0);
fprintf('Exact Vc                = %+.4f V\n',VcExact);
fprintf('Tlum(0 V)               = %.4f\n',Nom.Tlum0);
fprintf('DeltaTsol               = %.4f\n',Nom.DTsol);
fprintf('DeltaAsol(-10/+10 V)    = %+.4f / %+.4f\n',Nom.DAm,Nom.DAp);
fprintf('Reversal margin         = %.4f\n',Nom.Margin);
fprintf('====================================================\n');

%% One-at-a-time robustness
ParamName={'SiO_2^{(1)}','SiO_2^{(2)}','DAST','Au','Ag','EO coeff.'};
LevelName={'Low','Nominal','High'};
Np=6;Nl=3;
VcO=nan(Np,Nl);MarginO=zeros(Np,Nl);TlumO=zeros(Np,Nl);DTsolO=zeros(Np,Nl);
LowHigh=[85 90 95;110 115 120;95 100 105;4 5 6;7 8 9;.9 1 1.1];

fprintf('\nCalculating one-at-a-time robustness...\n');
for p=1:Np
    for l=1:Nl
        X=P0;
        switch p
            case 1,X(4)=LowHigh(p,l);
            case 2,X(5)=LowHigh(p,l);
            case 3,X(6:8)=LowHigh(p,l);
            case 4,X(1)=LowHigh(p,l);
            case 5,X(2:3)=LowHigh(p,l);
            case 6,X(9)=LowHigh(p,l);
        end
        R=evaluateDesign(X,wave,M,Vscan,ID,nBase,kDAST,dITO,dn_du,S,Veye,iSol,iLum);
        VcO(p,l)=R.Vc;MarginO(p,l)=R.Margin;TlumO(p,l)=R.Tlum0;DTsolO(p,l)=R.DTsol;
    end
    fprintf('%s completed\n',ParamName{p});
end

%% Monte Carlo: independent bounded fabrication/material variations
% Columns: Au, Ag1, Ag2, Si1, Si2, DAST1, DAST2, DAST3, EO scale
PMC=zeros(Nmc,9);
PMC(:,1)=randRange(4,6,Nmc);PMC(:,2)=randRange(7,9,Nmc);PMC(:,3)=randRange(7,9,Nmc);
PMC(:,4)=randRange(85,95,Nmc);PMC(:,5)=randRange(110,120,Nmc);
PMC(:,6)=randRange(95,105,Nmc);PMC(:,7)=randRange(95,105,Nmc);PMC(:,8)=randRange(95,105,Nmc);
PMC(:,9)=randRange(.9,1.1,Nmc);

VcMC=nan(Nmc,1);MarginMC=zeros(Nmc,1);TlumMC=zeros(Nmc,1);DTsolMC=zeros(Nmc,1);
DAmMC=zeros(Nmc,1);DApMC=zeros(Nmc,1);RecMC=zeros(Nmc,1);

fprintf('\nCalculating %d Monte Carlo structures...\n',Nmc);
for n=1:Nmc
    R=evaluateDesign(PMC(n,:),wave,M,Vscan,ID,nBase,kDAST,dITO,dn_du,S,Veye,iSol,iLum);
    VcMC(n)=R.Vc;MarginMC(n)=R.Margin;TlumMC(n)=R.Tlum0;DTsolMC(n)=R.DTsol;
    DAmMC(n)=R.DAm;DApMC(n)=R.DAp;RecMC(n)=R.Rec;
    if mod(n,25)==0||n==Nmc,fprintf('%4d / %d completed\n',n,Nmc);end
end

%% Robustness criteria
Reversal=DAmMC>0 & DApMC<0 & isfinite(VcMC);
Optical=TlumMC>=0.55 & DTsolMC>=0.018;
LowVoltage=abs(VcMC)<=7;
StrongMargin=MarginMC>=0.005;
FullPass=Reversal & Optical & LowVoltage & StrongMargin;

revRate=100*mean(Reversal);optRate=100*mean(Optical);passRate=100*mean(FullPass);
vr=VcMC(Reversal);mr=MarginMC(Reversal);

fprintf('\n================ FINAL ROBUSTNESS SUMMARY =================\n');
fprintf('Monte Carlo samples             = %d\n',Nmc);
fprintf('Directional-reversal persistence= %.1f %% (%d/%d)\n',revRate,sum(Reversal),Nmc);
fprintf('Smart-window optical criteria   = %.1f %%\n',optRate);
fprintf('Full robustness pass rate       = %.1f %% (%d/%d)\n',passRate,sum(FullPass),Nmc);
fprintf('Max reciprocity error           = %.3e\n',max(RecMC));

if ~isempty(vr)
    qVc=localPercentile(vr,[5 50 95]);
    qM=localPercentile(mr,[5 50 95]);
    fprintf('Vc 5/50/95 percentile           = %+.3f / %+.3f / %+.3f V\n',qVc);
    fprintf('Mrev 5/50/95 percentile         = %.4f / %.4f / %.4f\n',qM);
end
fprintf('Tlum(0 V) range                 = %.4f - %.4f\n',min(TlumMC),max(TlumMC));
fprintf('DeltaTsol range                 = %.4f - %.4f\n',min(DTsolMC),max(DTsolMC));
fprintf('DeltaAsol(-10 V) range          = %+.4f - %+.4f\n',min(DAmMC),max(DAmMC));
fprintf('DeltaAsol(+10 V) range          = %+.4f - %+.4f\n',min(DApMC),max(DApMC));

[~,iworst]=min(MarginMC);
fprintf('\nSmallest reversal-margin sample:\n');
fprintf('Au=%.2f, Ag1=%.2f, Ag2=%.2f, Si1=%.2f, Si2=%.2f nm\n',...
    PMC(iworst,1),PMC(iworst,2),PMC(iworst,3),PMC(iworst,4),PMC(iworst,5));
fprintf('DAST1/2/3=%.2f/%.2f/%.2f nm, EO scale=%.3f\n',...
    PMC(iworst,6),PMC(iworst,7),PMC(iworst,8),PMC(iworst,9));
fprintf('Vc=%+.3f V, DA(-10)=%+.4f, DA(+10)=%+.4f, Mrev=%.4f\n',...
    VcMC(iworst),DAmMC(iworst),DApMC(iworst),MarginMC(iworst));
fprintf('===========================================================\n');

%% Figure
figure('Color','w','Position',[30 50 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');ax=gobjects(1,8);

% (a) OAT Vc
ax(1)=nexttile;
imagesc(ax(1),1:3,1:Np,VcO);set(ax(1),'Color',[.85 .85 .85]);caxis(ax(1),[-10 10]);
set(ax(1),'XTick',1:3,'XTickLabel',LevelName,'YTick',1:Np,'YTickLabel',ParamName,'TickLabelInterpreter','tex');
cb=colorbar(ax(1));cb.Label.String='Reversal voltage {\it V}_c (V)';cb.Label.Interpreter='tex';

% (b) OAT reversal margin
ax(2)=nexttile;
imagesc(ax(2),1:3,1:Np,MarginO);
set(ax(2),'XTick',1:3,'XTickLabel',LevelName,'YTick',1:Np,'YTickLabel',ParamName,'TickLabelInterpreter','tex');
cb=colorbar(ax(2));cb.Label.String='Reversal margin {\it M}_{rev}';cb.Label.Interpreter='tex';

% (c) OAT Tlum0
ax(3)=nexttile;
imagesc(ax(3),1:3,1:Np,TlumO);caxis(ax(3),[0 1]);
set(ax(3),'XTick',1:3,'XTickLabel',LevelName,'YTick',1:Np,'YTickLabel',ParamName,'TickLabelInterpreter','tex');
cb=colorbar(ax(3));cb.Label.String='{\it T}_{lum}(0 V)';cb.Label.Interpreter='tex';

% (d) OAT DeltaTsol
ax(4)=nexttile;
imagesc(ax(4),1:3,1:Np,DTsolO);
set(ax(4),'XTick',1:3,'XTickLabel',LevelName,'YTick',1:Np,'YTickLabel',ParamName,'TickLabelInterpreter','tex');
cb=colorbar(ax(4));cb.Label.String='\Delta {\it T}_{sol}';cb.Label.Interpreter='tex';

% (e) Monte Carlo Vc-margin trade-off
ax(5)=nexttile;hold(ax(5),'on');
ok=Reversal;
scatter(ax(5),VcMC(ok),MarginMC(ok),30,TlumMC(ok),'filled');
plot(ax(5),VcExact,Nom.Margin,'kp','MarkerSize',12,'LineWidth',1.8);
xline(ax(5),0,'k:','HandleVisibility','off');
xlabel(ax(5),'Reversal voltage {\it V}_c (V)','Interpreter','tex');
ylabel(ax(5),'Reversal margin {\it M}_{rev}','Interpreter','tex');
cb=colorbar(ax(5));cb.Label.String='{\it T}_{lum}(0 V)';cb.Label.Interpreter='tex';
legend(ax(5),{'Monte Carlo','Nominal'},'Location','best','Box','off');

% (f) End-state polarity
ax(6)=nexttile;hold(ax(6),'on');
scatter(ax(6),DAmMC,DApMC,30,DTsolMC,'filled');
plot(ax(6),Nom.DAm,Nom.DAp,'kp','MarkerSize',12,'LineWidth',1.8);
xline(ax(6),0,'k:','HandleVisibility','off');yline(ax(6),0,'k:','HandleVisibility','off');
xlabel(ax(6),'\Delta {\it A}_{sol}(-10 V)','Interpreter','tex');
ylabel(ax(6),'\Delta {\it A}_{sol}(+10 V)','Interpreter','tex');
cb=colorbar(ax(6));cb.Label.String='\Delta {\it T}_{sol}';cb.Label.Interpreter='tex';
legend(ax(6),{'Monte Carlo','Nominal'},'Location','best','Box','off');

% (g) Vc distribution
ax(7)=nexttile;hold(ax(7),'on');
if any(Reversal),histogram(ax(7),VcMC(Reversal),18,'Normalization','probability');end
xline(ax(7),VcExact,'k--','LineWidth',LW);
xline(ax(7),0,'k:','HandleVisibility','off');
xlabel(ax(7),'Reversal voltage {\it V}_c (V)','Interpreter','tex');
ylabel(ax(7),'Probability');

% (h) Pass rates
ax(8)=nexttile;
Rates=[mean(Reversal) mean(Optical) mean(LowVoltage&Reversal) mean(StrongMargin&Reversal) mean(FullPass)]*100;
bar(ax(8),1:5,Rates,.6);
set(ax(8),'XTick',1:5,'XTickLabel',{'Reversal','Optical','|V_c|\leq7','M_{rev}\geq0.005','Full pass'},...
    'TickLabelInterpreter','tex');
xtickangle(ax(8),25);ylabel(ax(8),'Fraction (%)');ylim(ax(8),[0 105]);
for k=1:5,text(ax(8),k,Rates(k)+2,sprintf('%.1f%%',Rates(k)),'HorizontalAlignment','center','FontSize',FS1-3);end

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
function R=evaluateDesign(P,wave,M,Vscan,ID,nBase,kDAST,dITO,dn_du,S,Veye,iSol,iLum)
% P=[Au Ag1 Ag2 Si1 Si2 DAST1 DAST2 DAST3 EOscale]
dAu=P(1);dAg1=P(2);dAg2=P(3);dSi1=P(4);dSi2=P(5);tD=P(6:8);sc=P(9);
d_nm=[dAu tD(1) dITO dSi1 dITO tD(2) dAg1 dSi2 dAg2 tD(3) dITO];

Nv=numel(Vscan);DA=zeros(Nv,1);Tsol=zeros(Nv,1);R.Rec=0;R.Tlum0=NaN;
for iv=1:Nv
    V=Vscan(iv);

    [~,Tf,~,~,Tb,~]=calcState(wave,M,[V V V],d_nm,tD,ID,nBase,kDAST,dn_du,sc);
    Tsol(iv)=weighted(Tf,wave,S,iSol);
    R.Rec=max(R.Rec,max(abs(Tf-Tb)));
    if abs(V)<1e-12,R.Tlum0=weighted(Tf,wave,Veye,iLum);end

    [~,Tf1,Af1,~,Tb1,Ab1]=calcState(wave,M,[V 0 0],d_nm,tD,ID,nBase,kDAST,dn_du,sc);
    DA(iv)=weighted(Af1-Ab1,wave,S,iSol);
    R.Rec=max(R.Rec,max(abs(Tf1-Tb1)));
end
R.DA=DA;R.DAm=DA(1);R.DAp=DA(end);
R.Margin=min(R.DAm,-R.DAp);
R.Vc=zeroCrossNearest(Vscan,DA,0);
R.DTsol=max(Tsol)-min(Tsol);
end

function [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,V,d_nm,tD,ID,nBase,kDAST,dn_du,scale)
N=numel(wave);Rf=zeros(N,1);Tf=Rf;Af=Rf;Rb=Rf;Tb=Rf;Ab=Rf;df=d_nm*1e-9;
for iw=1:N
    d1=(nBase(iw)+scale*dn_du*V(1)/tD(1))-1i*kDAST(iw);
    d2=(nBase(iw)+scale*dn_du*V(2)/tD(2))-1i*kDAST(iw);
    d3=(nBase(iw)+scale*dn_du*V(3)/tD(3))-1i*kDAST(iw);
    nf=[M(iw,ID.Au),d1,M(iw,ID.ITO),M(iw,ID.SiO2),M(iw,ID.ITO),...
        d2,M(iw,ID.Ag),M(iw,ID.SiO2),M(iw,ID.Ag),d3,M(iw,ID.ITO)];
    lam=wave(iw)*1e-9;
    [Rf(iw),Tf(iw)]=tmmRT(1,1,nf,df,lam);
    [Rb(iw),Tb(iw)]=tmmRT(1,1,fliplr(nf),fliplr(df),lam);
    Af(iw)=1-Rf(iw)-Tf(iw);Ab(iw)=1-Rb(iw)-Tb(iw);
end
end

function val=solarDirectional(V1,P,wave,M,ID,nBase,kDAST,dITO,dn_du,S,iSol)
d_nm=[P(1) P(6) dITO P(4) dITO P(7) P(2) P(5) P(3) P(8) dITO];
[~,~,Af,~,~,Ab]=calcState(wave,M,[V1 0 0],d_nm,P(6:8),ID,nBase,kDAST,dn_du,P(9));
val=weighted(Af-Ab,wave,S,iSol);
end

function Vc=refineVc(fun,V,y)
Vc=zeroCrossNearest(V,y,0);
if ~isfinite(Vc),return;end
for k=1:numel(V)-1
    if y(k)*y(k+1)<0
        r=fzero(fun,[V(k) V(k+1)]);
        if abs(r)<=abs(Vc),Vc=r;end
    end
end
end

function Vc=zeroCrossNearest(V,y,Vref)
roots=[];
for k=1:numel(V)-1
    if abs(y(k))<1e-12,roots(end+1)=V(k);end
    if y(k)*y(k+1)<0
        roots(end+1)=V(k)-y(k)*(V(k+1)-V(k))/(y(k+1)-y(k));
    end
end
if abs(y(end))<1e-12,roots(end+1)=V(end);end
if isempty(roots),Vc=NaN;else,[~,i]=min(abs(roots-Vref));Vc=roots(i);end
end

function [R,T]=tmmRT(n0,ns,nf,df,lam)
MM=eye(2);
for j=1:numel(nf)
    n=nf(j);d=2*pi*n*df(j)/lam;
    MM=MM*[cos(d) 1i*sin(d)/n;1i*n*sin(d) cos(d)];
end
den=n0*MM(1,1)+n0*ns*MM(1,2)+MM(2,1)+ns*MM(2,2);
r=(n0*MM(1,1)+n0*ns*MM(1,2)-MM(2,1)-ns*MM(2,2))/den;
t=2*n0/den;
R=real(abs(r)^2);T=real(max(real(ns),0)/max(real(n0),eps)*abs(t)^2);
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

function y=weighted(x,w,W,idx)
x=x(:);w=w(:);W=W(:);
y=trapz(w(idx),x(idx).*W(idx))/trapz(w(idx),W(idx));
end

function x=randRange(a,b,N)
x=a+(b-a)*rand(N,1);
end

function q=localPercentile(x,p)
x=sort(x(:));N=numel(x);q=zeros(size(p));
for k=1:numel(p)
    r=1+(N-1)*p(k)/100;i=floor(r);f=r-i;
    if i>=N,q(k)=x(end);else,q(k)=x(i)+f*(x(i+1)-x(i));end
end
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