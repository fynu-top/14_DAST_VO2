% exportgraphics(gcf,'SiO2_spacer_optimization.pdf','ContentType','vector');
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

wave=(400:5:2500)';M=nk_data(wave);                 % n-i*k
ID.SiO2=2;ID.ITO=4;ID.DAST=10;ID.Au=11;ID.Ag=13;
assert(size(M,2)>=13,'nk_data must contain Ag as material ID 13.');
assert(all(isfinite(M(:,[ID.SiO2 ID.ITO ID.DAST ID.Au ID.Ag])),'all'),...
    'Material optical constants contain NaN/Inf over 400-2500 nm.');

%% Fixed structure
dAu=5;dDAST=100;dITO=40;dAg=8;                     % nm
n0=2.2;dn_du=3.41;                                  % nm/V
dSi1Scan=20:10:200;dSi2Scan=20:10:200;
Vopt=[-10 0 10];Vscan=(-10:1:10)';
baseSi=[80 90];MinBright=0.60;

N1=numel(dSi1Scan);N2=numel(dSi2Scan);Ndesign=N1*N2;

%% Spectral weights
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
assert(trapz(wave(iLum),Veye(iLum))>0&&trapz(wave(iSol),S(iSol))>0,'Invalid spectral weights.');

%% Coarse spacer scan
par=zeros(Ndesign,2);
Tbright=zeros(Ndesign,1);Tdark=Tbright;Tzero=Tbright;
DTsol=Tbright;DirCommon=Tbright;V1Span=Tbright;V1Reversal=false(Ndesign,1);RecErr=Tbright;
MapDT=zeros(N2,N1);MapDir=MapDT;MapT=MapDT;MapV1=MapDT;MapRev=false(N2,N1);

fprintf('\nScanning %d SiO2 spacer combinations...\n',Ndesign);
idx=0;
for i1=1:N1
    for i2=1:N2
        idx=idx+1;s1=dSi1Scan(i1);s2=dSi2Scan(i2);
        d_nm=[dAu dDAST dITO s1 dITO dDAST dAg s2 dAg dDAST dITO];

        tl=zeros(1,3);ts=tl;da=tl;da1=tl;rr=tl;
        for iv=1:3
            V=Vopt(iv);

            % Common voltage: V1=V2=V3
            [~,Tf,Af,~,Tb,Ab]=calcState(wave,M,[V V V],d_nm,dDAST,ID,n0,dn_du);
            tl(iv)=weighted(Tf,wave,Veye,iLum);ts(iv)=weighted(Tf,wave,S,iSol);
            da(iv)=weighted(Af,wave,S,iSol)-weighted(Ab,wave,S,iSol);
            rr(iv)=max(abs(Tf-Tb));

            % V1 only
            [~,~,Af1,~,~,Ab1]=calcState(wave,M,[V 0 0],d_nm,dDAST,ID,n0,dn_du);
            da1(iv)=weighted(Af1,wave,S,iSol)-weighted(Ab1,wave,S,iSol);
        end

        par(idx,:)=[s1 s2];
        Tbright(idx)=max(tl);Tdark(idx)=min(tl);Tzero(idx)=tl(2);
        DTsol(idx)=max(ts)-min(ts);DirCommon(idx)=max(abs(da));
        V1Span(idx)=max(da1)-min(da1);
        V1Reversal(idx)=min(da1)<0&&max(da1)>0;
        RecErr(idx)=max(rr);

        MapDT(i2,i1)=DTsol(idx);MapDir(i2,i1)=DirCommon(idx);
        MapT(i2,i1)=Tbright(idx);MapV1(i2,i1)=V1Span(idx);MapRev(i2,i1)=V1Reversal(idx);
    end
    fprintf('SiO2_1 = %3g nm completed\n',dSi1Scan(i1));
end

%% Feasible + Pareto search
feasible=Tbright>=MinBright;
F=find(feasible);
if isempty(F),warning('No structure satisfies Tlum,max >= %.2f; using all designs.',MinBright);F=(1:Ndesign)';end

PF=pareto2(F,DTsol,DirCommon);

% Balanced candidate: Pareto point closest to ideal normalized (1,1)
x=DTsol(PF);y=DirCommon(PF);
xn=(x-min(x))/max(max(x)-min(x),eps);
yn=(y-min(y))/max(max(y)-min(y),eps);
[~,kk]=min((1-xn).^2+(1-yn).^2);
iBalanced=PF(kk);

% Strongest V1 directional-reversal candidate
RF=F(V1Reversal(F));
if isempty(RF)
    fprintf('No V1 sign reversal found in the coarse {-10,0,+10 V} scan.\n');
    [~,kk]=max(V1Span(F));iReverse=F(kk);
else
    bestSpan=max(V1Span(RF));
    C=RF(abs(V1Span(RF)-bestSpan)<1e-12);
    [~,kk]=max(DTsol(C));iReverse=C(kk);
end

ibase=find(par(:,1)==baseSi(1)&par(:,2)==baseSi(2),1);
if isempty(ibase),error('Baseline [80 90] nm is not contained in the scan.');end

fprintf('\n================ SPACER OPTIMIZATION =================\n');
fprintf('Total designs    : %d\n',Ndesign);
fprintf('Tlum-feasible    : %d\n',numel(F));
fprintf('Pareto designs   : %d\n',numel(PF));
fprintf('V1 reversal      : %d designs\n',sum(V1Reversal));
fprintf('Reciprocity error: %.3e\n',max(RecErr));

printCandidate('Baseline',ibase,par,Tbright,Tdark,Tzero,DTsol,DirCommon,V1Span,V1Reversal);
printCandidate('Balanced Pareto',iBalanced,par,Tbright,Tdark,Tzero,DTsol,DirCommon,V1Span,V1Reversal);
printCandidate('V1 reversal',iReverse,par,Tbright,Tdark,Tzero,DTsol,DirCommon,V1Span,V1Reversal);
fprintf('======================================================\n');

%% Full voltage verification
idxDesign=[ibase iBalanced iReverse];
Name={'Baseline','Balanced','V1 reversal'};
Nd=3;Nv=numel(Vscan);

Tlum=zeros(Nv,Nd);Tsol=zeros(Nv,Nd);DAcommon=zeros(Nv,Nd);DAV1=zeros(Nv,Nd);
TfShow=cell(Nd,3);Vshow=[-10 0 10];

for d=1:Nd
    s1=par(idxDesign(d),1);s2=par(idxDesign(d),2);
    d_nm=[dAu dDAST dITO s1 dITO dDAST dAg s2 dAg dDAST dITO];

    for iv=1:Nv
        V=Vscan(iv);

        [~,Tf,Af,~,~,Ab]=calcState(wave,M,[V V V],d_nm,dDAST,ID,n0,dn_du);
        Tlum(iv,d)=weighted(Tf,wave,Veye,iLum);
        Tsol(iv,d)=weighted(Tf,wave,S,iSol);
        DAcommon(iv,d)=weighted(Af,wave,S,iSol)-weighted(Ab,wave,S,iSol);

        [~,~,Af1,~,~,Ab1]=calcState(wave,M,[V 0 0],d_nm,dDAST,ID,n0,dn_du);
        DAV1(iv,d)=weighted(Af1,wave,S,iSol)-weighted(Ab1,wave,S,iSol);

        q=find(Vshow==V,1);if ~isempty(q),TfShow{d,q}=Tf;end
    end
end

fprintf('\n================ FULL VOLTAGE CHECK ==================\n');
for d=1:Nd
    s1=par(idxDesign(d),1);s2=par(idxDesign(d),2);
    Vcross=zeroCross(Vscan,DAV1(:,d));
    fprintf('%s: SiO2_1=%g nm, SiO2_2=%g nm\n',Name{d},s1,s2);
    fprintf('  Tlum max/min    = %.4f / %.4f\n',max(Tlum(:,d)),min(Tlum(:,d)));
    fprintf('  Delta Tsol      = %.4f\n',max(Tsol(:,d))-min(Tsol(:,d)));
    fprintf('  max|DeltaAsol|  = %.4f\n',max(abs(DAcommon(:,d))));
    fprintf('  V1 span Asol    = %.4f\n',max(DAV1(:,d))-min(DAV1(:,d)));
    if isfinite(Vcross),fprintf('  V1 reversal near V = %.2f V\n',Vcross);
    else,fprintf('  V1 reversal     = not found over -10:+10 V\n');end
end
fprintf('======================================================\n');

%% Figure
figure('Color','w','Position',[30 60 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');ax=gobjects(1,8);

% (a) Delta Tsol map
ax(1)=nexttile;
imagesc(ax(1),dSi1Scan,dSi2Scan,MapDT);axis(ax(1),'xy');
xlabel(ax(1),'SiO_2 spacer 1 (nm)');ylabel(ax(1),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(1));cb.Label.String='\Delta {\it T}_{sol}';cb.Label.Interpreter='tex';

% (b) Common-voltage directional solar contrast
ax(2)=nexttile;
imagesc(ax(2),dSi1Scan,dSi2Scan,MapDir);axis(ax(2),'xy');
xlabel(ax(2),'SiO_2 spacer 1 (nm)');ylabel(ax(2),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(2));cb.Label.String='max|\Delta {\it A}_{sol}|';cb.Label.Interpreter='tex';

% (c) Maximum luminous transmittance
ax(3)=nexttile;
imagesc(ax(3),dSi1Scan,dSi2Scan,MapT);axis(ax(3),'xy');caxis(ax(3),[0 1]);
xlabel(ax(3),'SiO_2 spacer 1 (nm)');ylabel(ax(3),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(3));cb.Label.String='max {\it T}_{lum}';cb.Label.Interpreter='tex';

% (d) V1 directional-control span
ax(4)=nexttile;
imagesc(ax(4),dSi1Scan,dSi2Scan,MapV1);axis(ax(4),'xy');
xlabel(ax(4),'SiO_2 spacer 1 (nm)');ylabel(ax(4),'SiO_2 spacer 2 (nm)');
cb=colorbar(ax(4));cb.Label.String='V_1 span of \Delta {\it A}_{sol}';cb.Label.Interpreter='tex';
hold(ax(4),'on');
[r,c]=find(MapRev);
if ~isempty(r),plot(ax(4),dSi1Scan(c),dSi2Scan(r),'ko','MarkerSize',5,'LineWidth',1);end

% (e) Pareto plane
ax(5)=nexttile;hold(ax(5),'on');
sc=scatter(ax(5),DTsol(F),DirCommon(F),36,Tbright(F),'filled');
scatter(ax(5),DTsol(PF),DirCommon(PF),65,'o','LineWidth',1.2);
scatter(ax(5),DTsol(ibase),DirCommon(ibase),100,'s','LineWidth',1.8);
scatter(ax(5),DTsol(iBalanced),DirCommon(iBalanced),130,'p','LineWidth',2);
scatter(ax(5),DTsol(iReverse),DirCommon(iReverse),130,'d','LineWidth',2);
xlabel(ax(5),'\Delta {\it T}_{sol}','Interpreter','tex');
ylabel(ax(5),'max|\Delta {\it A}_{sol}|','Interpreter','tex');
cb=colorbar(ax(5));cb.Label.String='max {\it T}_{lum}';cb.Label.Interpreter='tex';
legend(ax(5),{'Feasible','Pareto','Baseline','Balanced','V_1 reversal'},...
    'Location','best','Box','off','Interpreter','tex','FontSize',FS1-4);

% (f) Full common-voltage Tlum
ax(6)=nexttile;hold(ax(6),'on');LS={'--','-','-.'};
for d=1:Nd,plot(ax(6),Vscan,Tlum(:,d),LS{d},'LineWidth',LW);end
xlabel(ax(6),'Common voltage (V)');ylabel(ax(6),'{\it T}_{lum}','Interpreter','tex');
legend(ax(6),Name,'Location','best','Box','off');xlim(ax(6),[-10 10]);ylim(ax(6),[0 1]);

% (g) Full common-voltage directional response
ax(7)=nexttile;hold(ax(7),'on');
for d=1:Nd,plot(ax(7),Vscan,DAcommon(:,d),LS{d},'LineWidth',LW);end
yline(ax(7),0,'k:','HandleVisibility','off');
xlabel(ax(7),'Common voltage (V)');ylabel(ax(7),'\Delta {\it A}_{sol}','Interpreter','tex');
legend(ax(7),Name,'Location','best','Box','off');xlim(ax(7),[-10 10]);

% (h) V1-only directional reversal
ax(8)=nexttile;hold(ax(8),'on');
for d=1:Nd,plot(ax(8),Vscan,DAV1(:,d),LS{d},'LineWidth',LW);end
yline(ax(8),0,'k:','HandleVisibility','off');
xlabel(ax(8),'{\it V}_1 only (V)','Interpreter','tex');
ylabel(ax(8),'\Delta {\it A}_{sol}','Interpreter','tex');
legend(ax(8),Name,'Location','best','Box','off');xlim(ax(8),[-10 10]);

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
    Af(iw)=1-Rf(iw)-Tf(iw);Ab(iw)=1-Rb(iw)-Tb(iw);
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
r=(n0*M(1,1)+n0*ns*M(1,2)-M(2,1)-ns*M(2,2))/den;t=2*n0/den;
R=real(abs(r)^2);T=real(max(real(ns),0)/max(real(n0),eps)*abs(t)^2);
end

function y=weighted(x,w,W,idx)
x=x(:);w=w(:);W=W(:);
y=trapz(w(idx),x(idx).*W(idx))/trapz(w(idx),W(idx));
end

function PF=pareto2(F,x,y)
keep=true(size(F));
for a=1:numel(F)
    i=F(a);
    for b=1:numel(F)
        j=F(b);
        if j~=i&&x(j)>=x(i)&&y(j)>=y(i)&&(x(j)>x(i)||y(j)>y(i))
            keep(a)=false;break;
        end
    end
end
PF=F(keep);
end

function printCandidate(name,i,par,Tbright,Tdark,Tzero,DTsol,Dir,V1Span,Rev)
fprintf('\n%s: SiO2_1=%g nm, SiO2_2=%g nm\n',name,par(i,1),par(i,2));
fprintf('  Tlum,max       = %.4f\n',Tbright(i));
fprintf('  Tlum,min       = %.4f\n',Tdark(i));
fprintf('  Tlum(0 V)      = %.4f\n',Tzero(i));
fprintf('  Delta Tsol     = %.4f\n',DTsol(i));
fprintf('  max|DeltaAsol| = %.4f\n',Dir(i));
fprintf('  V1 Asol span   = %.4f\n',V1Span(i));
fprintf('  V1 reversal    = %s\n',string(Rev(i)));
end

function Vc=zeroCross(V,y)
Vc=NaN;
for k=1:numel(V)-1
    if y(k)==0,Vc=V(k);return;end
    if y(k)*y(k+1)<0
        Vc=V(k)-y(k)*(V(k+1)-V(k))/(y(k+1)-y(k));return;
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
        f=fullfile(folders{k},names{j});if isfile(f),file=f;return;end
    end
    f=which(names{j});if ~isempty(f),file=f;return;end
end
error('Cannot find data file: %s',strjoin(names,', '));
end

function [w,y]=readWeight(file,col)
M=readmatrix(file);assert(size(M,2)>=col,'Insufficient columns in %s.',file);
ok=isfinite(M(:,1))&isfinite(M(:,col))&M(:,1)>0;w=M(ok,1);y=M(ok,col);
if max(w)<100,w=w*1000;end
[w,id]=sort(w);y=y(id);[w,id]=unique(w,'stable');y=y(id);
end