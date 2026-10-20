
function Y = obj1(x,shoulian,nk_M,wave,IncidentAngle)

h = x;

h = [0 h 0*100]; % [ 上层PMMA 中间层 底层基底];
%%
clear Absorption xz
% wave = load('wave.mat');wave = wave.wave';
wave = wave';
%% case 4
%{ 1
R1 = zeros(length(wave),1); T1 = R1; A1 = R1;
for m = 1:length(wave)
    n0 = 1;
    nsubs = 1; % 基底材质为空气
    n_layers = nk_M(m,:);
    d_layers = h*1e-9;
    theta_in = 0;
%       [R1(m,1) T1(m,1) A1(m,1)] = tm1( wave(m), nk_M(m,:),  h  ,IncidentAngle);    
%       [R2(m,1) T2(m,1) A2(m,1)] = tm1( wave(m), rot90(rot90(nk_M(m,:))),  rot90(rot90(h))  ,IncidentAngle);

%     [z11 z12 z13] = RTA( wave(m), nk_M(m,:),  h  ,IncidentAngle, 0); % polarization
%     [z21 z22 z23] = RTA( wave(m), nk_M(m,:),  h  ,IncidentAngle, 1); % polarization
%     R1(m,1) = (z11 + z21)/2;    T1(m,1) = (z12 + z22)/2;    A1(m,1) = (z13 + z23)/2;    % TE+TM
    
    [R1(m,1) T1(m,1) A1(m,1)] = RTA( wave(m), nk_M(m,:),  h  ,IncidentAngle, 0); % TE and TM
%     [R2(m,1) T2(m,1) A2(m,1)] = RTA( wave(m), rot90(rot90(nk_M(m,:))),  rot90(rot90(h))  ,IncidentAngle,0);
end
    R2 = R1; T2 = T1; A2 = A1;
% figure('color',[1 1 1]); 
% plot(wave_range,R1,'b',  wave_range,T1,'r',  wave_range,A1,'g',  'LineWidth',2);hold on;
% plot(wave_range,R2,'b--',wave_range,T2,'r--',wave_range,A2,'g--','LineWidth',2);hold on;
% legend('\fontsize{20}R','\fontsize{20}T','\fontsize{20}A'); xlabel('wavelength [nm]');ylabel('intensity');
%% calculation
% V = interp1(Visible(:,1),Visible(:,2),wave)';
% Visible_I = V*T1/sum(V);
% Target=[540 660 640 960];
% for ix = 1:length(Target)
%     [ipx(ix) ipy(ix)] = find(wave == Target(ix));
% end
% Visible_intensity = sum(T1(ipx(1):ipx(2)))/length(T1(ipx(1):ipx(2)));
% A_difference =abs(sum(A1(ipx(3):ipx(4)))/length(A1(ipx(3):ipx(4)))-sum(A2(ipx(3):ipx(4)))/length(A2(ipx(3):ipx(4))) );

% if shoulian == 1     Y = -(  Visible_intensity+ A_difference*3     )  ; end
% if shoulian == 2     
Y = [R1 R2 T1 T2 A1 A2];  %end

% objectivefunction([61 32 9.2 46.3])[R1 T1 A1 R2 T2 A2]


