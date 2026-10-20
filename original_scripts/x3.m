% exportgraphics(gcf,'DAST_thickness_final_scan.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

FS1=20;FS2=26;FN='Times New Roman';LW=1.5;FigWidth=1900;FigHeight=880;
set(groot,'defaultAxesFontName',FN,'defaultTextFontName',FN,'defaultLegendFontName',FN,...
    'defaultAxesFontSize',FS1,'defaultTextFontSize',FS1,'defaultLegendFontSize',FS1,...
    'defaultAxesFontWeight','normal','defaultTextFontWeight','normal',...
    'defaultAxesLabelFontSizeMultiplier',1,'defaultAxesTitleFontSizeMultiplier',1,...
    'defaultTextInterpreter','tex','defaultLegendInterpreter','tex',...
    'defaultAxesTickLabelInterpreter','tex');

%% Paths and material data
root=fileparts(mfilename('fullpath'));if isempty(root),root=pwd;end
project=fileparts(root);
P={root,project,fullfile(root,'optical_constants'),fullfile(root,'script'),...
    fullfile(root,'Atmospheric data'),fullfile(project,'optical_constants'),...
    fullfile(project,'script'),fullfile(project,'Atmospheric data')};
for k=1:numel(P),if isfolder(P{k}),addpath(P{k});end,end

wave=(400:5:2500)';M=nk_data(wave);                 % complex index: n-i*k
ID.SiO2=2;ID.ITO=4;ID.DAST=10;ID.Au=11;ID.Ag=13;
assert(size(M,2)>=13,'nk_data must contain material IDs through Ag = 13.');
assert(all(isfinite(M(:,[ID.SiO2 ID.ITO ID.DAST ID.Au ID.Ag])),'all'),...
    'Material optical constants contain NaN/Inf over 400-2500 nm.');

%% Fixed structure
% Air/Au/DAST1/ITO/SiO2_1/ITO/DAST2/Ag/SiO2_2/Ag/DAST3/ITO/Air
dAu=5;dITO=40;dAg=8;dSi1=80;dSi2=90;                % nm
n0=2.2;dn_du=3.41;                                   % nm/V
tScan=(70:5:150)';Vscan=(-10:1:10)';
Nt=numel(tScan);Nv=numel(Vscan);

%% Spectral weights
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
assert(trapz(wave(iLum),Veye(iLum))>0&&trapz(wave(iSol),S(iSol))>0,'Invalid spectral weights.');

%% Allocate
Tlum=zeros(Nv,Nt);Tsol=Tlum;AsolF=Tlum;AsolB=Tlum;DAcommon=Tlum;
DAV1=Tlum;RecErr=Tlum;

%% Full thickness-voltage scan
fprintf('\nCalculating uniform DAST thickness scan...\n');
for it=1:Nt
    t=tScan(it);
    d_nm=[dAu t dITO dSi1 dITO t dAg dSi2 dAg t dITO];

    for iv=1:Nv
        V=Vscan(iv);

        % Common voltage: V1=V2=V3
        [~,Tf,Af,~,Tb,Ab]=calcState(wave,M,[V V V],d_nm,t,ID,n0,dn_du);
        Tlum(iv,it)=weighted(Tf,wave,Veye,iLum);
        Tsol(iv,it)=weighted(Tf,wave,S,iSol);
        AsolF(iv,it)=weighted(Af,wave,S,iSol);
        AsolB(iv,it)=weighted(Ab,wave,S,iSol);
        DAcommon(iv,it)=AsolF(iv,it)-AsolB(iv,it);
        RecErr(iv,it)=max(abs(Tf-Tb));

        % V1-only control
        [~,~,Af1,~,~,Ab1]=calcState(wave,M,[V 0 0],d_nm,t,ID,n0,dn_du);
        DAV1(iv,it)=weighted(Af1,wave,S,iSol)-weighted(Ab1,wave,S,iSol);
    end
    fprintf('tDAST = %3g nm completed\n',t);
end

%% Thickness-dependent metrics
[TlumMax,iMax]=max(Tlum,[],1);[TlumMin,iMin]=min(Tlum,[],1);
Tlum0=Tlum(Vscan==0,:);
TsolMax=max(Tsol,[],1);TsolMin=min(Tsol,[],1);
DeltaTlum=TlumMax-TlumMin;DeltaTsol=TsolMax-TsolMin;
MaxDir=max(abs(DAcommon),[],1);
V1Span=max(DAV1,[],1)-min(DAV1,[],1);

VmaxLum=Vscan(iMax);VminLum=Vscan(iMin);
Vc=nan(1,Nt);
for it=1:Nt,Vc(it)=zeroCross(Vscan,DAV1(:,it));end

DeltaN10=dn_du*10./tScan';                            % |Delta n| at 10 V
E10=10000./tScan';                                    % V/um at 10 V

%% Smart-window constraints
feasible=TlumMax>=0.60&TlumMin>=0.40&Tlum0>=0.55;
F=find(feasible);

[bestDT,iDT]=max(DeltaTsol);
[bestDir,iDir]=max(MaxDir);
[bestBright,iBright]=max(TlumMax);

fprintf('\n================ DAST THICKNESS SCREENING =================\n');
fprintf('Maximum reciprocity error = %.3e\n',max(RecErr(:)));
fprintf('Feasible thicknesses      = %d / %d\n',numel(F),Nt);

if ~isempty(F)
    fprintf('Feasible range            = ');
    fprintf('%g ',tScan(F));
    fprintf('nm\n');
else
    fprintf('No thickness satisfies all three smart-window constraints.\n');
end

fprintf('\nLargest DeltaTsol   = %.4f at t = %.0f nm\n',bestDT,tScan(iDT));
fprintf('Largest |DeltaAsol| = %.4f at t = %.0f nm\n',bestDir,tScan(iDir));
fprintf('Largest Tlum,max    = %.4f at t = %.0f nm\n',bestBright,tScan(iBright));

fprintf('\n tDAST  Tlum_max  Tlum_min  Tlum(0V)  DeltaTlum  DeltaTsol  max|DAsol|  Vc(V1)  |Dn|@10V  E@10V\n');
checkT=90:5:110;
for q=1:numel(checkT)
    it=find(tScan==checkT(q),1);
    if isfinite(Vc(it)),vcText=sprintf('%+.2f',Vc(it));else,vcText='  -- ';end
    fprintf('%5.0f   %.4f    %.4f     %.4f      %.4f      %.4f      %.4f     %6s    %.4f      %.1f\n',...
        tScan(it),TlumMax(it),TlumMin(it),Tlum0(it),DeltaTlum(it),...
        DeltaTsol(it),MaxDir(it),vcText,DeltaN10(it),E10(it));
end

i100=find(tScan==100,1);
fprintf('\n100-nm reference:\n');
fprintf('  n(-10 V)=%.4f, n(0 V)=%.4f, n(+10 V)=%.4f\n',...
    n0-DeltaN10(i100),n0,n0+DeltaN10(i100));
fprintf('  Tlum,max/min/0V = %.4f / %.4f / %.4f\n',TlumMax(i100),TlumMin(i100),Tlum0(i100));
fprintf('  DeltaTsol       = %.4f\n',DeltaTsol(i100));
fprintf('  max|DeltaAsol|  = %.4f\n',MaxDir(i100));
if isfinite(Vc(i100)),fprintf('  V1 reversal     = %.2f V\n',Vc(i100));
else,fprintf('  V1 reversal     = not found\n');end
fprintf('===========================================================\n');

%% Figure
figure('Color','w','Position',[30 60 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');ax=gobjects(1,8);

% (a) Tlum map
ax(1)=nexttile;
imagesc(ax(1),tScan,Vscan,Tlum);axis(ax(1),'xy');caxis(ax(1),[0 1]);
xlabel(ax(1),'DAST thickness (nm)');ylabel(ax(1),'Common voltage (V)');
cb=colorbar(ax(1));cb.Label.String='{\it T}_{lum}';cb.Label.Interpreter='tex';

% (b) Tsol map
ax(2)=nexttile;
imagesc(ax(2),tScan,Vscan,Tsol);axis(ax(2),'xy');caxis(ax(2),[0 1]);
xlabel(ax(2),'DAST thickness (nm)');ylabel(ax(2),'Common voltage (V)');
cb=colorbar(ax(2));cb.Label.String='{\it T}_{sol}';cb.Label.Interpreter='tex';

% (c) Common-voltage directional solar absorptance
ax(3)=nexttile;
m=max(abs(DAcommon),[],'all');if m<1e-6,m=1e-6;end
imagesc(ax(3),tScan,Vscan,DAcommon);axis(ax(3),'xy');caxis(ax(3),[-m m]);
xlabel(ax(3),'DAST thickness (nm)');ylabel(ax(3),'Common voltage (V)');
cb=colorbar(ax(3));cb.Label.String='\Delta {\it A}_{sol}';cb.Label.Interpreter='tex';

% (d) V1-only directional solar absorptance
ax(4)=nexttile;
m=max(abs(DAV1),[],'all');if m<1e-6,m=1e-6;end
imagesc(ax(4),tScan,Vscan,DAV1);axis(ax(4),'xy');caxis(ax(4),[-m m]);
xlabel(ax(4),'DAST thickness (nm)');ylabel(ax(4),'{\it V}_1 only (V)','Interpreter','tex');
cb=colorbar(ax(4));cb.Label.String='\Delta {\it A}_{sol}';cb.Label.Interpreter='tex';

% (e) Luminous-transmittance metrics
ax(5)=nexttile;hold(ax(5),'on');
plot(ax(5),tScan,TlumMax,'-o','LineWidth',LW,'MarkerSize',5);
plot(ax(5),tScan,Tlum0,'-s','LineWidth',LW,'MarkerSize',5);
plot(ax(5),tScan,TlumMin,'-^','LineWidth',LW,'MarkerSize',5);
xline(ax(5),100,':','LineWidth',1.2,'HandleVisibility','off');
xlabel(ax(5),'DAST thickness (nm)');ylabel(ax(5),'{\it T}_{lum}','Interpreter','tex');
legend(ax(5),{'Maximum','{\it V}=0 V','Minimum'},'Location','best','Box','off','Interpreter','tex');
xlim(ax(5),[tScan(1) tScan(end)]);ylim(ax(5),[0 1]);

% (f) Solar modulation + directional contrast
ax(6)=nexttile;hold(ax(6),'on');
yyaxis(ax(6),'left');
plot(ax(6),tScan,DeltaTsol,'-o','LineWidth',LW,'MarkerSize',5);
ylabel(ax(6),'\Delta {\it T}_{sol}','Interpreter','tex');
yyaxis(ax(6),'right');
plot(ax(6),tScan,MaxDir,'-s','LineWidth',LW,'MarkerSize',5);
ylabel(ax(6),'max|\Delta {\it A}_{sol}|','Interpreter','tex');
xline(ax(6),100,':','LineWidth',1.2,'HandleVisibility','off');
xlabel(ax(6),'DAST thickness (nm)');
legend(ax(6),{'\Delta {\it T}_{sol}','max|\Delta {\it A}_{sol}|'},...
    'Location','best','Box','off','Interpreter','tex');
xlim(ax(6),[tScan(1) tScan(end)]);

% (g) V1 reversal voltage and control span
ax(7)=nexttile;hold(ax(7),'on');
yyaxis(ax(7),'left');
plot(ax(7),tScan,Vc,'-o','LineWidth',LW,'MarkerSize',5);
yline(ax(7),0,':','LineWidth',1,'HandleVisibility','off');
ylabel(ax(7),'Reversal voltage {\it V}_c (V)','Interpreter','tex');
yyaxis(ax(7),'right');
plot(ax(7),tScan,V1Span,'-s','LineWidth',LW,'MarkerSize',5);
ylabel(ax(7),'{\it V}_1 span of \Delta {\it A}_{sol}','Interpreter','tex');
xline(ax(7),100,':','LineWidth',1.2,'HandleVisibility','off');
xlabel(ax(7),'DAST thickness (nm)');
xlim(ax(7),[tScan(1) tScan(end)]);

% (h) Electro-optic strength
ax(8)=nexttile;hold(ax(8),'on');
yyaxis(ax(8),'left');
plot(ax(8),tScan,DeltaN10,'-o','LineWidth',LW,'MarkerSize',5);
ylabel(ax(8),'max|\Delta {\it n}_{DAST}|','Interpreter','tex');
yyaxis(ax(8),'right');
plot(ax(8),tScan,E10,'-s','LineWidth',LW,'MarkerSize',5);
ylabel(ax(8),'Maximum electric field (V/\mum)','Interpreter','tex');
xline(ax(8),100,':','LineWidth',1.2,'HandleVisibility','off');
xlabel(ax(8),'DAST thickness (nm)');
legend(ax(8),{'|\Delta {\it n}| at 10 V','|{\it E}| at 10 V'},...
    'Location','best','Box','off','Interpreter','tex');
xlim(ax(8),[tScan(1) tScan(end)]);

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