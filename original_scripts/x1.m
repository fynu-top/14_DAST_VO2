% exportgraphics(gcf,'DAST_Janus_direct_electrode_test.pdf','ContentType','vector');
clear;clc;close all;addpath('optical_constants','script');tic;

FS1=20;FS2=26;FN='Times New Roman';LW=1.5;FigWidth=1900;FigHeight=860;mu=char(181);
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

wave=(400:5:2500)';M=nk_data(wave);                 % complex index: n-i*k
ID.SiO2=2;ID.ITO=4;ID.DAST=10;ID.Au=11;ID.Ag=13;
assert(size(M,2)>=13,'nk_data must contain Ag as material ID 13.');
assert(all(isfinite(M(:,[ID.SiO2 ID.ITO ID.DAST ID.Au ID.Ag])),'all'),...
    'Material optical constants contain NaN/Inf over 400-2500 nm.');

%% Structure
% Air/Au/DAST1/ITO/SiO2/ITO/DAST2/Ag/SiO2/Ag/DAST3/ITO/Air
dAu=5;dDAST=[100 100 100];dITO=[40 40 40];dSi=[80 90];dAg=[8 8]; % nm
d_nm=[dAu dDAST(1) dITO(1) dSi(1) dITO(2) dDAST(2) dAg(1) dSi(2) dAg(2) dDAST(3) dITO(3)];

%% DAST electro-optic model
n0=2.2;dn_du=3.41;                                  % nm/V
Vscan=(-10:1:10)';Vshow=[-10 0 10];

fprintf('\n================ DAST ELECTRODE TOPOLOGY ================\n');
fprintf('DAST1: Au / DAST1 / ITO\n');
fprintf('DAST2: ITO / DAST2 / Ag\n');
fprintf('DAST3: Ag / DAST3 / ITO\n');
for j=1:3
    fprintf('DAST%d: t=%g nm, n(-10V)=%.4f, n(0)=%.4f, n(+10V)=%.4f\n',...
        j,dDAST(j),n0-10*dn_du/dDAST(j),n0,n0+10*dn_du/dDAST(j));
end
fprintf('==========================================================\n');

%% Spectral weights
fileSol=findData(root,{'ASTMG173.xlsx','AM0AM1_5.xlsx','AM0AM1_5.csv'});
fileLum=findData(root,{'CIE_Vlambda.csv','CIE_Vlambda.xlsx'});
[ws,S0]=readWeight(fileSol,3);[wv,V0]=readWeight(fileLum,2);
S=interp1(ws,S0,wave,'linear',0);Veye=interp1(wv,V0,wave,'linear',0);
iLum=wave>=400&wave<=780;iSol=wave>=400&wave<=2500;
assert(trapz(wave(iLum),Veye(iLum))>0&&trapz(wave(iSol),S(iSol))>0,'Invalid spectral weights.');

%% Common-voltage sweep
Nv=numel(Vscan);Nw=numel(wave);
RfMap=zeros(Nv,Nw);TfMap=RfMap;AfMap=RfMap;RbMap=RfMap;TbMap=RfMap;AbMap=RfMap;
Tlum=zeros(Nv,1);Tsol=Tlum;AsolF=Tlum;AsolB=Tlum;RsolF=Tlum;RsolB=Tlum;RecErr=Tlum;

fprintf('\nCalculating common-voltage sweep...\n');
for iv=1:Nv
    VV=[Vscan(iv) Vscan(iv) Vscan(iv)];
    [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,VV,d_nm,dDAST,ID,n0,dn_du);
    RfMap(iv,:)=Rf;TfMap(iv,:)=Tf;AfMap(iv,:)=Af;
    RbMap(iv,:)=Rb;TbMap(iv,:)=Tb;AbMap(iv,:)=Ab;
    Tlum(iv)=weighted(Tf,wave,Veye,iLum);Tsol(iv)=weighted(Tf,wave,S,iSol);
    AsolF(iv)=weighted(Af,wave,S,iSol);AsolB(iv)=weighted(Ab,wave,S,iSol);
    RsolF(iv)=weighted(Rf,wave,S,iSol);RsolB(iv)=weighted(Rb,wave,S,iSol);
    RecErr(iv)=max(abs(Tf-Tb));
end

%% Independent DAST control
TlumInd=zeros(Nv,3);TsolInd=zeros(Nv,3);AsolDirInd=zeros(Nv,3);
fprintf('Calculating independent V1/V2/V3 control...\n');
for iv=1:Nv
    for j=1:3
        VV=[0 0 0];VV(j)=Vscan(iv);
        [~,Tf,Af,~,~,Ab]=calcState(wave,M,VV,d_nm,dDAST,ID,n0,dn_du);
        TlumInd(iv,j)=weighted(Tf,wave,Veye,iLum);
        TsolInd(iv,j)=weighted(Tf,wave,S,iSol);
        AsolDirInd(iv,j)=weighted(Af,wave,S,iSol)-weighted(Ab,wave,S,iSol);
    end
end

%% Main metrics
DeltaAsol=AsolF-AsolB;
DeltaTlum=max(Tlum)-min(Tlum);
DeltaTsol=max(Tsol)-min(Tsol);
MaxDirSol=max(abs(DeltaAsol));
AdirMap=AfMap-AbMap;
MaxDirSpec=max(abs(AdirMap),[],'all');

fprintf('\n================ SMART-WINDOW RESULTS =================\n');
fprintf('Tlum max/min          = %.4f / %.4f\n',max(Tlum),min(Tlum));
fprintf('Delta Tlum            = %.4f (%.2f percentage points)\n',DeltaTlum,100*DeltaTlum);
fprintf('Tsol max/min          = %.4f / %.4f\n',max(Tsol),min(Tsol));
fprintf('Delta Tsol            = %.4f (%.2f percentage points)\n',DeltaTsol,100*DeltaTsol);
fprintf('max |Delta Asol|      = %.4f\n',MaxDirSol);
fprintf('max spectral |Af-Ab|  = %.4f\n',MaxDirSpec);
fprintf('max reciprocity error = %.3e\n',max(RecErr));

for Vq=Vshow
    iv=find(Vscan==Vq,1);
    fprintf('\nV = %+g V\n',Vq);
    fprintf('  Tlum       = %.4f\n',Tlum(iv));
    fprintf('  Tsol       = %.4f\n',Tsol(iv));
    fprintf('  Asol,F     = %.4f\n',AsolF(iv));
    fprintf('  Asol,B     = %.4f\n',AsolB(iv));
    fprintf('  Delta Asol = %+.4f\n',DeltaAsol(iv));
end

fprintf('\nIndependent channels, -10 V to +10 V:\n');
for j=1:3
    fprintf('V%d only: DeltaTlum=%.4f, DeltaTsol=%.4f, max|DeltaAsol|=%.4f\n',...
        j,max(TlumInd(:,j))-min(TlumInd(:,j)),...
        max(TsolInd(:,j))-min(TsolInd(:,j)),max(abs(AsolDirInd(:,j))));
end
fprintf('=======================================================\n');

%% Figure
figure('Color','w','Position',[30 60 FigWidth FigHeight]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');ax=gobjects(1,8);C=lines(3);

% (a) Transmission spectra
ax(1)=nexttile;hold(ax(1),'on');
for q=1:3
    iv=find(Vscan==Vshow(q),1);
    plot(ax(1),wave/1000,TfMap(iv,:),'LineWidth',LW,'Color',C(q,:));
end
xlabel(ax(1),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(1),'{\it T}','Interpreter','tex');
legend(ax(1),{'-10 V','0 V','+10 V'},'Location','best','Box','off');
xlim(ax(1),[.4 2.5]);ylim(ax(1),[0 1]);

% (b) Forward/reverse absorptance
ax(2)=nexttile;hold(ax(2),'on');
for q=1:3
    iv=find(Vscan==Vshow(q),1);
    plot(ax(2),wave/1000,AfMap(iv,:),'-','LineWidth',LW,'Color',C(q,:));
    plot(ax(2),wave/1000,AbMap(iv,:),'--','LineWidth',LW,'Color',C(q,:));
end
xlabel(ax(2),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(2),'{\it A}','Interpreter','tex');
legend(ax(2),{'F,-10 V','B,-10 V','F,0 V','B,0 V','F,+10 V','B,+10 V'},...
    'Location','best','NumColumns',2,'Box','off','FontSize',FS1-4);
xlim(ax(2),[.4 2.5]);ylim(ax(2),[0 1]);

% (c) Voltage-transmission map
ax(3)=nexttile;
imagesc(ax(3),wave/1000,Vscan,TfMap);axis(ax(3),'xy');caxis(ax(3),[0 1]);
xlabel(ax(3),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(3),'Common voltage (V)');
cb=colorbar(ax(3));cb.Label.String='{\it T}';cb.Label.Interpreter='tex';

% (d) Directional absorptance map
ax(4)=nexttile;
m=max(abs(AdirMap),[],'all');if m<1e-6,m=1e-6;end
imagesc(ax(4),wave/1000,Vscan,AdirMap);axis(ax(4),'xy');caxis(ax(4),[-m m]);
xlabel(ax(4),['Wavelength (' mu 'm)'],'Interpreter','none');
ylabel(ax(4),'Common voltage (V)');
cb=colorbar(ax(4));cb.Label.String='\Delta {\it A}_{dir}={\it A}_{F}-{\it A}_{B}';cb.Label.Interpreter='tex';

% (e) Weighted transmittance
ax(5)=nexttile;hold(ax(5),'on');
plot(ax(5),Vscan,Tlum,'-o','LineWidth',LW,'MarkerSize',5);
plot(ax(5),Vscan,Tsol,'-s','LineWidth',LW,'MarkerSize',5);
xlabel(ax(5),'Common voltage (V)');
ylabel(ax(5),'Weighted transmittance');
legend(ax(5),{'{\it T}_{lum}','{\it T}_{sol}'},'Location','best','Box','off','Interpreter','tex');
xlim(ax(5),[-10 10]);ylim(ax(5),[0 1]);

% (f) Forward/reverse solar absorption
ax(6)=nexttile;hold(ax(6),'on');
yyaxis(ax(6),'left');
plot(ax(6),Vscan,AsolF,'-o','LineWidth',LW,'MarkerSize',5);
plot(ax(6),Vscan,AsolB,'-s','LineWidth',LW,'MarkerSize',5);
ylabel(ax(6),'{\it A}_{sol}','Interpreter','tex');
yyaxis(ax(6),'right');
plot(ax(6),Vscan,DeltaAsol,'-^','LineWidth',LW,'MarkerSize',5);
ylabel(ax(6),'\Delta {\it A}_{sol}','Interpreter','tex');
xlabel(ax(6),'Common voltage (V)');
legend(ax(6),{'Forward','Reverse','Difference'},'Location','best','Box','off');
xlim(ax(6),[-10 10]);

% (g) Independent Tsol channels
ax(7)=nexttile;hold(ax(7),'on');
plot(ax(7),Vscan,TsolInd(:,1),'-','LineWidth',LW);
plot(ax(7),Vscan,TsolInd(:,2),'--','LineWidth',LW);
plot(ax(7),Vscan,TsolInd(:,3),'-.','LineWidth',LW);
plot(ax(7),Vscan,Tsol,':','LineWidth',LW+0.3);
xlabel(ax(7),'Applied voltage (V)');
ylabel(ax(7),'{\it T}_{sol}','Interpreter','tex');
legend(ax(7),{'{\it V}_1 only','{\it V}_2 only','{\it V}_3 only','{\it V}_1={\it V}_2={\it V}_3'},...
    'Location','best','Box','off','Interpreter','tex');
xlim(ax(7),[-10 10]);ylim(ax(7),[0 1]);

% (h) Independent directional response
ax(8)=nexttile;hold(ax(8),'on');
plot(ax(8),Vscan,AsolDirInd(:,1),'-','LineWidth',LW);
plot(ax(8),Vscan,AsolDirInd(:,2),'--','LineWidth',LW);
plot(ax(8),Vscan,AsolDirInd(:,3),'-.','LineWidth',LW);
plot(ax(8),Vscan,DeltaAsol,':','LineWidth',LW+0.3);
yline(ax(8),0,'k:','HandleVisibility','off');
xlabel(ax(8),'Applied voltage (V)');
ylabel(ax(8),'\Delta {\it A}_{sol}','Interpreter','tex');
legend(ax(8),{'{\it V}_1 only','{\it V}_2 only','{\it V}_3 only','Common voltage'},...
    'Location','best','Box','off','Interpreter','tex');
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
function [Rf,Tf,Af,Rb,Tb,Ab]=calcState(wave,M,V,d_nm,dDAST,ID,n0,dn_du)
N=numel(wave);Rf=zeros(N,1);Tf=Rf;Af=Rf;Rb=Rf;Tb=Rf;Ab=Rf;
for iw=1:N
    d1=dastIndex(M(iw,ID.DAST),V(1),dDAST(1),n0,dn_du);
    d2=dastIndex(M(iw,ID.DAST),V(2),dDAST(2),n0,dn_du);
    d3=dastIndex(M(iw,ID.DAST),V(3),dDAST(3),n0,dn_du);
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

function file=findData(root,names)
project=fileparts(root);
folders={root,project,fullfile(root,'Atmospheric data'),fullfile(project,'Atmospheric data')};
file='';
for j=1:numel(names)
    for k=1:numel(folders)
        f=fullfile(folders{k},names{j});
        if isfile(f),file=f;return;end
    end
    f=which(names{j});if ~isempty(f),file=f;return;end
end
error('Cannot find data file: %s',strjoin(names,', '));
end

function [w,y]=readWeight(file,col)
M=readmatrix(file);assert(size(M,2)>=col,'Insufficient columns in %s.',file);
ok=isfinite(M(:,1))&isfinite(M(:,col))&M(:,1)>0;
w=M(ok,1);y=M(ok,col);
if max(w)<100,w=w*1000;end
[w,id]=sort(w);y=y(id);[w,id]=unique(w,'stable');y=y(id);
end