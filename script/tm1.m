% load('aTMM.mat')
function [R T A] = tm1(lambda, nk,  dM ,IncidentAngle)
nk=conj(nk);
nk = [1 nk 1];
dM = [0 dM 0];

thetaM = zeros(1,length(nk));
deltaM = thetaM;
etaM = thetaM;
Matrices = zeros(2,2,length(nk));

% nk = [ n_enter  nAg  n_silica nTDBC n_silica nAg nS];
% dM =     [0       30     50      20      50    30  0];
thetaM(1) = IncidentAngle;  % incident angle
deltaM(1) = 2*pi*nk(1)*dM(1)*cos(thetaM(1))/lambda;
etaM(1) = nk(1)*cos(thetaM(1));
Matrices = [cos(deltaM(1)),-1i*sin(deltaM(1))/etaM(1); % primera matriz
           -1i*etaM(1)*sin(deltaM(1)),cos(deltaM(1))];
Mcon = Matrices;
%{ 1
for im = 2:length(nk)
    thetaM(im) = 0;%asin( nk(im-1)/nk(im)   *sin(thetaM(im-1))   );
    deltaM(im) = 2*pi*nk(im)*dM(im)*cos(thetaM(im))/lambda;
    etaM(im)   = nk(im)*cos(thetaM(im));
    Matrices(:,:,im) = [cos(deltaM(im)),-1i*sin(deltaM(im))/etaM(im); % matrices de cada capa
                       -1i*etaM(im)*sin(deltaM(im)),cos(deltaM(im))];
	Mcon = Mcon*Matrices(:,:,im);
end
%}
%%
%{ 
    thetaM1 = thetaM;
    deltaM1 = 2*pi*nk.*dM.*cos(thetaM1)/lambda;
    etaM1   = nk.*cos(thetaM1);
    Matrices1(:,:,1:length(nk)) = [cos(deltaM),-1i*sin(deltaM)./etaM; % matrices de cada capa
                       -1i*etaM.*sin(deltaM),cos(deltaM)];
                   
    
    im = 2:length(nk);
    %thetaM(im) = asin( nk(im-1)/nk(im)   *sin(thetaM(im-1))   );
    thetaM(1:length(nk)) = 0;
    deltaM(im) = 2*pi*nk(im).*dM(im).*cos(thetaM(im))./lambda;
    etaM(im)   = nk(im)*cos(thetaM(im));
    Matrices(:,:,im) = [cos(deltaM(im)),-1i*sin(deltaM(im))/etaM(im); % matrices de cada capa
                       -1i*etaM(im)*sin(deltaM(im)),cos(deltaM(im))];
%     Mcon = Mcon*Matrices(:,:,im);
    Mcon = prod(Mcon);
%}
%%
etaS = nk(im)*cos(thetaM(im));
eta_m = etaS;
eta_o = nk(1)*cos(thetaM(1));
BC = Mcon*[1;etaS];
Y = BC(2)/BC(1);
R = ((eta_o-Y)/(eta_o+Y))*conj((eta_o-Y)/(eta_o+Y)); % reflectancia
A = (4*eta_o*real(BC(1)*conj(BC(2))-eta_m))/((eta_o*BC(1)+BC(2))*conj(eta_o*BC(1)+BC(2))); % absorbancia
T = (4*eta_o*real(eta_m))/((eta_o*BC(1)+BC(2))*conj(eta_o*BC(1)+BC(2))); % transmitancia
% if T==0     A = 1-R; end

end

