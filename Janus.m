classdef Janus
% Shared implementation of the uploaded n-i*k TMM and DAST EO model.
methods (Static)
function D=init(step,physicalSiO2,photopicCSV)
if nargin<1,step=5;end;if nargin<2,physicalSiO2=false;end;if nargin<3,photopicCSV=false;end
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'script'),fullfile(root,'optical_constants'));
D.root=root;D.w=(400:step:2500)';D.M=nk_data(D.w);
if physicalSiO2
 X=readmatrix(fullfile(root,'optical_constants','2SiO2.xlsx'));D.M(:,2)=interp1(X(:,1)*1000,X(:,2)-1i*X(:,3),D.w);
end
assert(all(isfinite(D.M(:,[2 4 10 11 13])),'all'),'Material data contain NaN or Inf.');
D.n=Janus.base(D.w);D.k=abs(imag(D.M(:,10)));D.coeff=3.41;
D.d=[5 100 40 90 40 100 8 115 8 100 40];D.ids=[11 10 4 2 4 10 13 2 13 10 4];
D.labels={'Au','DAST_1','ITO_1','SiO_2^{(1)}','ITO_2','DAST_2','Ag_1','SiO_2^{(2)}','Ag_2','DAST_3','ITO_3'};
A=readmatrix(fullfile(root,'VO2_temperature_data','ASTMG173.xlsx'));ok=isfinite(A(:,1))&isfinite(A(:,3))&A(:,1)>0;[w,id]=unique(A(ok,1));v=A(ok,3);D.S=interp1(w,v(id),D.w,'linear',0);
% CIE 1931 observer: [wavelength,xbar,ybar,zbar] -> V(lambda)=ybar (column 3).
% A genuine two-column V(lambda) file is also accepted when photopicCSV=true.
if photopicCSV,cieFile=fullfile(root,'VO2_temperature_data','CIE_Vlambda.csv');else,cieFile=fullfile(root,'Atmospheric data','CIE_Vlambda.xlsx');end
assert(isfile(cieFile),'Missing CIE photopic file: %s',cieFile);A=readmatrix(cieFile);
if size(A,2)>=4,col=3;elseif size(A,2)==2,col=2;else,error('CIE file must have [w,V] or [w,xbar,ybar,zbar] columns.');end
ok=isfinite(A(:,1))&isfinite(A(:,col))&A(:,1)>=300& A(:,1)<=900;
assert(nnz(ok)>100,'Insufficient CIE photopic data.');[w,id]=unique(A(ok,1));a=A(ok,col);D.Vphotopic=interp1(w,a(id),D.w,'linear',0);
% Obtain the official CIE standard illuminant D65 CSV (wavelength_nm, relative_SPD).
% https://www.cie.co.at/datatable/cie-standard-illuminant-d65
D65file=fullfile(root,'Atmospheric data','CIE_std_illum_D65.csv');
if ~isfile(D65file),D65file=fullfile(root,'Atmospheric data','CIE_std_illum_D65.xlsx');end
assert(isfile(D65file),'Missing official D65 data: place CIE_std_illum_D65.csv in Atmospheric data/ (see CIE DOI 10.25039/CIE.DS.hjfjmt59).');
B=readmatrix(D65file);ok=isfinite(B(:,1))&isfinite(B(:,2))&B(:,1)>=300&B(:,1)<=900;
assert(nnz(ok)>100,'Invalid/insufficient D65 spectrum.');[w,id]=unique(B(ok,1));b=B(ok,2);D.D65=interp1(w,b(id),D.w,'linear',0);
D.eye=D.Vphotopic.*D.D65;D.Wlum=D.eye;D.lum=D.w>=400&D.w<=780;
D.sol=true(size(D.w));D.physicalSiO2=physicalSiO2;D.photopicCSV=photopicCSV;
assert(all(isfinite(D.eye(D.lum)))&&all(D.eye(D.lum)>=0)&&trapz(D.w(D.lum),D.eye(D.lum))>0,'Invalid D65*ybar weights.');
assert(trapz(D.w,D.S)>0,'Invalid solar spectral weights.');
end
function n=base(w)
n=2.2*ones(size(w));i=w>650;ns=sqrt(2.078^2+1.645*533^2./(w(i).^2-533^2));s=ones(size(ns));j=w(i)<750;wi=w(i);s(j)=.5-.5*cos(pi*(wi(j)-650)/100);n(i)=(1-s)*2.2+s.*ns;
end
function N=stack(D,V)
N=D.M(:,D.ids);for j=1:3,k=[2 6 10];N(:,k(j))=D.n+D.coeff*V(j)/D.d(k(j))-1i*D.k;end
end
function kz=kz(n,s)
kz=sqrt(n.^2-s.^2);i=real(kz)<0|(abs(real(kz))<1e-12&imag(kz)>0);kz(i)=-kz(i);
end
function [R,T,r,t]=rt(N,d,w,theta,pol)
if nargin<4,theta=0;end;if nargin<5,pol='TE';end
if theta>=90,R=ones(size(w));T=zeros(size(w));r=-ones(size(w));t=zeros(size(w));return;end
s=sind(theta);eta0=cosd(theta);if strcmpi(pol,'TM'),eta0=1/eta0;end
A=ones(size(w));B=zeros(size(w));C=B;E=A;
for j=1:numel(d)
 n=N(:,j);kz=Janus.kz(n,s);eta=kz;if strcmpi(pol,'TM'),eta=n.^2./kz;end
 delta=2*pi*kz*d(j)./w;c=cos(delta);u=1i*sin(delta)./eta;v=1i*eta.*sin(delta);
 a=A.*c+B.*v;b=A.*u+B.*c;c2=C.*c+E.*v;e=C.*u+E.*c;A=a;B=b;C=c2;E=e;
end
bb=A+B*eta0;cc=C+E*eta0;den=eta0*bb+cc;r=(eta0*bb-cc)./den;t=2*eta0./den;R=abs(r).^2;T=abs(t).^2;
end
function Q=state(D,V,theta,pol)
if nargin<3,theta=0;end;if nargin<4,pol='TE';end
N=Janus.stack(D,V);[Q.RF,Q.TF]=Janus.rt(N,D.d,D.w,theta,pol);[Q.RB,Q.TB]=Janus.rt(fliplr(N),fliplr(D.d),D.w,theta,pol);Q.AF=1-Q.RF-Q.TF;Q.AB=1-Q.RB-Q.TB;Q.DA=Q.AF-Q.AB;
Q.Tlum=Janus.weight(D,Q.TF,'lum');Q.Tsol=Janus.weight(D,Q.TF,'sol');Q.AsF=Janus.weight(D,Q.AF,'sol');Q.AsB=Janus.weight(D,Q.AB,'sol');Q.DAsol=Q.AsF-Q.AsB;Q.rec=max(abs(Q.TF-Q.TB));
end
function y=weight(D,x,type)
if strcmp(type,'lum'),i=D.lum;W=D.eye;else,i=D.sol;W=D.S;end
y=trapz(D.w(i),x(i,:).*W(i),1)/trapz(D.w(i),W(i));
end
function R=sweep(D,V,ch,theta,pol)
if nargin<4,theta=0;end;if nargin<5,pol='TE';end
R.V=V(:);R.Tlum=zeros(numel(V),1);R.Tsol=R.Tlum;R.AsF=R.Tlum;R.AsB=R.Tlum;R.DA=R.Tlum;R.rec=0;
for j=1:numel(V),u=zeros(1,3);if ch==0,u(:)=V(j);else,u(ch)=V(j);end;Q=Janus.state(D,u,theta,pol);R.Tlum(j)=Q.Tlum;R.Tsol(j)=Q.Tsol;R.AsF(j)=Q.AsF;R.AsB(j)=Q.AsB;R.DA(j)=Q.DAsol;R.rec=max(R.rec,Q.rec);end
end
function r=roots(V,y)
r=[];for k=1:numel(V)-1,if abs(y(k))<1e-10,r(end+1)=V(k);end;if y(k)*y(k+1)<0,r(end+1)=V(k)-y(k)*(V(k+1)-V(k))/(y(k+1)-y(k));end;end;if abs(y(end))<1e-10,r(end+1)=V(end);end;r=unique(round(r,8));
end
function c=vc(D)
V=(-10:1:10)';R=Janus.sweep(D,V,1);r=Janus.roots(V,R.DA);if isempty(r),c=NaN;return;end;[~,j]=min(abs(r));c=r(j);k=find(V<=c,1,'last');if k<numel(V)&&R.DA(k)*R.DA(k+1)<0,c=fzero(@(v)Janus.da(D,v),[V(k) V(k+1)]);end
end
function y=da(D,v)
Q=Janus.state(D,[v 0 0]);y=Q.DAsol;
end
function R=design(D,V)
C=Janus.sweep(D,V,0);S=Janus.sweep(D,V,1);r=Janus.roots(V,S.DA);R.Vc=NaN;if ~isempty(r),[~,j]=min(abs(r));R.Vc=r(j);end
Q=Janus.state(D,[0 0 0]);R.Tlum0=Q.Tlum;R.DTsol=max(C.Tsol)-min(C.Tsol);R.maxDA=max(abs(S.DA));R.DAm=S.DA(1);R.DAp=S.DA(end);R.margin=min(R.DAm,-R.DAp);R.rec=max(C.rec,S.rec);
end
function [R,T,a,v]=layerRT(n,d,w)
N=numel(n);J=cell(N,1);M=eye(2);for k=1:N,delta=2*pi*n(k)*d(k)/w;J{k}=[cos(delta) 1i*sin(delta)/n(k);1i*n(k)*sin(delta) cos(delta)];M=M*J{k};end
b=M(1,1)+M(1,2);c=M(2,1)+M(2,2);r=(b-c)/(b+c);t=2/(b+c);R=abs(r)^2;T=abs(t)^2;v=[1+r;1-r];a=zeros(1,N);
for k=1:N,p=real(v(1)*conj(v(2)));u=J{k}\v;a(k)=p-real(u(1)*conj(u(2)));v=u;end
a(abs(a)<1e-12)=0;
end
function Q=layers(D,V)
Q=Janus.state(D,V);N=Janus.stack(D,V);Q.LF=zeros(numel(D.w),numel(D.d));Q.LB=Q.LF;
for j=1:numel(D.w),[~,~,Q.LF(j,:)]=Janus.layerRT(N(j,:),D.d,D.w(j));[~,~,a]=Janus.layerRT(fliplr(N(j,:)),fliplr(D.d),D.w(j));Q.LB(j,:)=fliplr(a);end
Q.LsolF=Janus.weight(D,Q.LF,'sol');Q.LsolB=Janus.weight(D,Q.LB,'sol');Q.LDA=Q.LsolF-Q.LsolB;
end
function [z,F,B]=field(D,V,lambda,np)
if nargin<4,np=40;end;N=Janus.stack(D,V);n=interp1(D.w,N,lambda,'linear');[z,F]=Janus.oneField(n,D.d,lambda,np);[zb,eb]=Janus.oneField(fliplr(n),fliplr(D.d),lambda,np);[zb,i]=sort(sum(D.d)-zb);B=interp1(zb,eb(i),z,'linear');
end
function [z,e]=oneField(n,d,w,np)
M=eye(2);for k=1:numel(n),delta=2*pi*n(k)*d(k)/w;M=M*[cos(delta) 1i*sin(delta)/n(k);1i*n(k)*sin(delta) cos(delta)];end
b=M(1,1)+M(1,2);c=M(2,1)+M(2,2);r=(b-c)/(b+c);v=[1+r;1-r];z=[];e=[];z0=0;
for k=1:numel(n)
 x=linspace(0,d(k),np+1);if k>1,x=x(2:end);end
 for j=1:numel(x),delta=2*pi*n(k)*x(j)/w;u=[cos(delta) 1i*sin(delta)/n(k);1i*n(k)*sin(delta) cos(delta)]\v;z(end+1,1)=z0+x(j);e(end+1,1)=abs(u(1))^2;end
 delta=2*pi*n(k)*d(k)/w;v=[cos(delta) 1i*sin(delta)/n(k);1i*n(k)*sin(delta) cos(delta)]\v;z0=z0+d(k);
end
end
function [fig,ax]=canvas(nr,nc,count,width,height)
fig=figure('Color','w','Position',[30 40 width height]);tiledlayout(fig,nr,nc,'TileSpacing','loose','Padding','loose');ax=gobjects(1,count);for j=1:count,ax(j)=nexttile;hold(ax(j),'on');end
end
function lg=leg(ax,labels,varargin)
lg=legend(ax,labels,'Location','best','Box','off','AutoUpdate','off',varargin{:});
end
function finish(fig,ax,FS,FS2,FN,LW,gaps)
drawnow;for j=1:numel(ax)
 set(ax(j),'FontName',FN,'FontSize',FS,'FontWeight','normal','LineWidth',LW,'TickDir','in','TickLength',[.012 .012],'Box','on','Layer','top');set([ax(j).XLabel ax(j).YLabel],'FontWeight','normal','FontName',FN,'FontSize',FS);
 text(ax(j),0,1.025,sprintf('(%s)',char(96+j)),'Units','normalized','FontName',FN,'FontSize',FS2,'FontWeight','normal','HorizontalAlignment','left','VerticalAlignment','bottom','Clipping','off','HandleVisibility','off');
end
lg=gobjects(0,1);for j=1:numel(ax),if ~isempty(ax(j).Legend)&&isvalid(ax(j).Legend),lg(end+1,1)=ax(j).Legend;end;end;assert(size(gaps,1)==numel(lg),'Provide one xGap/yGap pair for each legend.');drawnow;
for j=1:numel(lg),set(lg(j),'FontName',FN,'FontSize',FS,'FontWeight','normal','Units','normalized');p=lg(j).Position;lg(j).Position=p+[gaps(j,:) 0 0];end
cb=findall(fig,'Type','colorbar');for j=1:numel(cb),set(cb(j),'FontName',FN,'FontSize',FS,'LineWidth',LW,'TickDirection','in');cb(j).Label.FontName=FN;cb(j).Label.FontSize=FS;cb(j).Label.FontWeight='normal';cb(j).Label.Interpreter='tex';end
end
function heat(ax,x,y,z,label,diverging)
if nargin<6,diverging=false;end
imagesc(ax,x,y,z);axis(ax,'xy');xlim(ax,[min(x) max(x)]);ylim(ax,[min(y) max(y)]);cb=colorbar(ax);cb.Label.String=label;
if diverging,colormap(ax,Janus.divmap());a=max(abs(z(:)),[],'omitnan');if isfinite(a)&&a>0,clim(ax,[-a a]);end;else,colormap(ax,parula);end
end
function C=divmap()
x=[0 .5 1];c=[.16 .32 .7;1 1 1;.75 .12 .16];C=interp1(x,c,linspace(0,1,256));
end
function name=states(vc)
name={'-10 V',sprintf('{\\it V}_c=%.3f V',vc),'+10 V'};
end
function a=continuous(theta,mask)
j=find(~mask,1);if isempty(j),a=theta(end);elseif j==1,a=NaN;else,a=theta(j-1);end
end
function D=fromP(D,p)
D.d=[p(1) p(6) 40 p(4) 40 p(7) p(2) p(5) p(3) p(8) 40];D.coeff=3.41*p(9);
end
function R=mc(D,N,V)
% Original independent uniform distributions and rng(1).
rng(1);lo=[4 7 7 85 110 95 95 95 .9];hi=[6 9 9 95 120 105 105 105 1.1];P=zeros(N,9);for j=1:9,P(:,j)=lo(j)+(hi(j)-lo(j))*rand(N,1);end
R.P=P;R.Vc=nan(N,1);R.margin=R.Vc;R.Tlum0=R.Vc;R.DTsol=R.Vc;R.DAm=R.Vc;R.DAp=R.Vc;R.rec=R.Vc;
for j=1:N,Q=Janus.design(Janus.fromP(D,P(j,:)),V);f={'Vc','margin','Tlum0','DTsol','DAm','DAp','rec'};for k=1:numel(f),R.(f{k})(j)=Q.(f{k});end;if mod(j,25)==0,fprintf('Monte Carlo: %d/%d\n',j,N);end;end
R.reversal=R.DAm>0&R.DAp<0&isfinite(R.Vc);R.optical=R.Tlum0>=.55&R.DTsol>=.018;R.lowV=abs(R.Vc)<=7;R.strong=R.margin>=.005;R.full=R.reversal&R.optical&R.lowV&R.strong;
end
function A=angles(D,theta,V)
A.theta=theta(:);A.V=V(:);A.Vc0=Janus.vc(D);A.vstate=[-10 A.Vc0 10];na=numel(theta);nv=numel(V);A.DA=zeros(na,nv,2);A.Tlum=A.DA;A.Tsol=A.DA;A.AsF=zeros(na,3,2);A.AsB=A.AsF;A.TlumState=A.AsF;A.TsolState=A.AsF;A.DAState=A.AsF;A.rec=0;
for p=1:2,pol={'TE','TM'};for j=1:na
 R=Janus.sweep(D,V,1,theta(j),pol{p});A.DA(j,:,p)=R.DA.';A.Tlum(j,:,p)=R.Tlum.';A.Tsol(j,:,p)=R.Tsol.';A.rec=max(A.rec,R.rec);
 for s=1:3,Q=Janus.state(D,[A.vstate(s) 0 0],theta(j),pol{p});A.AsF(j,s,p)=Q.AsF;A.AsB(j,s,p)=Q.AsB;A.TlumState(j,s,p)=Q.Tlum;A.TsolState(j,s,p)=Q.Tsol;A.DAState(j,s,p)=Q.DAsol;end
end;fprintf('%s angle map complete\n',pol{p});end
A.UDA=mean(A.DA,3);A.UTlum=mean(A.TlumState,3);A.UTsol=mean(A.TsolState,3);A.UDAState=mean(A.DAState,3);
end
function G=maps(D,x,y,kind,V)
G.x=x;G.y=y;names={'Vc','margin','DAm','DAp','Tlum0','DTsol','maxDA'};for j=1:numel(names),G.(names{j})=nan(numel(y),numel(x));end
for i=1:numel(x),for j=1:numel(y),E=D;if strcmp(kind,'spacer'),E.d(4)=x(i);E.d(8)=y(j);else,E.d(1)=x(i);E.d([7 9])=y(j);end;Q=Janus.design(E,V);for k=1:numel(names),G.(names{k})(j,i)=Q.(names{k});end;end;fprintf('Scan: %d/%d columns complete\n',i,numel(x));end
end
function O=oat(D,V)
O.names={'SiO_2^{(1)}','SiO_2^{(2)}','DAST','Au','Ag','EO coeff.'};levels=[85 90 95;110 115 120;95 100 105;4 5 6;7 8 9;.9 1 1.1];O.Vc=nan(6,3);O.margin=O.Vc;O.Tlum0=O.Vc;O.DTsol=O.Vc;
for j=1:6,for k=1:3,E=D;switch j,case 1,E.d(4)=levels(j,k);case 2,E.d(8)=levels(j,k);case 3,E.d([2 6 10])=levels(j,k);case 4,E.d(1)=levels(j,k);case 5,E.d([7 9])=levels(j,k);case 6,E.coeff=3.41*levels(j,k);end;Q=Janus.design(E,V);f={'Vc','margin','Tlum0','DTsol'};for m=1:4,O.(f{m})(j,k)=Q.(f{m});end;end;end
end

end
end
