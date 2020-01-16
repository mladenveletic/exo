step = 1;          %[ms]
npulse = 1;
time = 0:step:1e3; %[ms]
amp = 20;          %[\muA/cm^2; pA]
offset = 250;
%% Parameter Set: ---
% NEURONAL CELLS
% Nominal values for electrical activity
% Voltages given in [mV]
% Conductances given in [nS or, equivalently, mS/cm^3]
% Time constants given in [ms]
% Capacitance given in [\muF/cm^2]
gK = 36;
gNa = 120;
gL = 0.3;
VNa = 50;
VK = -70;
VL = -54.4;
cm = 1;

VCa = 65; 
VSOC = 0; 
gSOC = 0;
gCaL = 0.7; 
gCaPQ = 0.6;
gCaT = 0.4;
gKa = 1;
gKATP = 0.295;
VmCaL = -30;
VmCaPQ = -5; 
SmCaL = 10;
SmCaPQ = 10; 
VhCaL = -33; 
VhCaPQ = -33; 
ShCaL = -5;
ShCaPQ = -5;
taumVCaL = 1;
taumVCaPQ = 1;
taum0CaL = 0.05;
taum0CaPQ = 0.05;
VtaumCaL = -23;
VtaumCaPQ = -23;
StaumCaL = 20;
StaumCaPQ = 20;
tauhVCaL = 60;
tauhVCaPQ = 60;
tauh0CaL = 51;
tauh0CaPQ = 51;
VtauhCaL = 0;
VtauhCaPQ = 0;
StauhCaL = 20;
StauhCaPQ = 20;
VmCaT = -49;
SmCaT = 4;
VhCaT = -52;
ShCaT = -5;
taumVCaT = 15;
taum0CaT = 0;
VtaumCaT = -50;
StaumCaT = 12;
tauhVCaT = 20;
tauh0CaT = 5;
VtauhCaT = -50;
StauhCaT = 15;

VmNa = -30;
SmNa = 4;
VhNa = -52;
ShNa = -8;
taumVNa = 6;
taum0Na = 0.05;
VtaumNa = -50;
StaumNa = 10;
tauhVNa = 120;
tauh0Na = 0.5;
VtauhNa = -50;
StauhNa = 8;
VmKa = -45;
SmKa = 10;
VhKa = -68;
ShKa = -10;
taumVKa = 0;
taum0Ka = 0.1;
tauhVKa = 60;
tauh0Ka = 5;
VtauhKa = 5; 
StauhKa = 20;
VmK = -25;
SmK = 23;
taumVK = 1.5;
taum0K = 15;
VtaumK = -10;
StaumK = 25;

f = 0.01;
Volmud = 2.618e-19; % [L]
Volc = 5.725e-13;   % [L]
NPQ = 0.0001;
Bmud = 264;         % [1/ms]
kPMCA = 0.3;        % [1/ms]
pleak = 3e-4;       % [1/ms]
nPQ = 4;
nL = 4; 
nm = 4;
alpha = 5.18e-15;   % [\mumol/pA/ms]
Volm = 5.149e-14;   % [L]
VolcDIVVoler = 31;
NL = 200;           
Bm = 0.128;         % [1/ms]
kSERCA = 0.1;       % [1/ms]
KPQ = 2;
KL = 50;            % [muM]
Km = 2;             % [muM]
%% Derivatives
% x(1) = Vm;
% x(2) = mCaL;
% x(3) = mCaPQ;
% x(4) = mCaT;
% x(5) = mNa;   -->m
% x(6) = mK;    -->n
% x(7) = mKa;
% x(8) = hCaL;
% x(9) = hCaPQ;
% x(10) = hCaT;
% x(11) = hNa;  -->h
% x(12) = hKa;
% x(13) = CaL;
% x(14) = CaPQ;
% x(15) = Cam;
% x(16) = Cac;
% x(17) = Caer;
%-------------------------------------------------------------------------- 
f1 = @(t,x) [
-(gK*x(6)^4*(x(1)-VK) + gNa*x(5)^3*x(11)*(x(1)-VNa) + gL*(x(1)-VL)- controlSignal(t,amp,time,step,npulse, offset))/cm;
(funcminf(x(1), VmCaL, SmCaL) - x(2))./functaum(x(1), taumVCaL, VtaumCaL, StaumCaL, taum0CaL);
(funcminf(x(1), VmCaPQ, SmCaPQ) - x(3))./functaum(x(1), taumVCaPQ, VtaumCaPQ, StaumCaPQ, taum0CaPQ);
(funcminf(x(1), VmCaT, SmCaT) - x(4))./functaum(x(1), taumVCaT, VtaumCaT, StaumCaT, taum0CaT);
0.1*(x(1)+40)/(1-exp(-(x(1)+40)/10))*(1-x(5))-4*exp(-(x(1)+65)/18)*x(5);
0.01*(x(1)+55)/(1-exp(-(x(1)+55)/10))*(1-x(6))-0.125*exp(-(x(1)+65)/80)*x(6);
(funcminf(x(1), VmKa, SmKa) - x(7))./functaum(x(1), taumVKa, VtaumK, StaumK, taum0Ka);
(funcminf(x(1), VhCaL, ShCaL) - x(8))./functaum(x(1), tauhVCaL, VtauhCaL, StauhCaL, tauh0CaL);
(funcminf(x(1), VhCaPQ, ShCaPQ) - x(9))./functaum(x(1), tauhVCaPQ, VtauhCaPQ, StauhCaPQ, tauh0CaPQ);
(funcminf(x(1), VhCaT, ShCaT) - x(10))./functaum(x(1), tauhVCaT, VtauhCaT, StauhCaT, tauh0CaT);
0.07*exp(-(x(1)+65)/20)*(1-x(11))-1/(1+exp(-(x(1)+35)/10))*x(11);
(funcminf(x(1), VhKa, ShKa) - x(12))./functaum(x(1), tauhVKa, VtauhKa, StauhKa, tauh0Ka);
(-f*alpha*gCaL.*(x(1)-VCa)./NL./Volmud - f*Bmud*(x(13)-x(15)));
(-f*alpha*gCaPQ.*(x(1)-VCa)./NPQ./Volmud*0 - f*Bmud*(x(14)-x(15)));
(f*(-alpha*(gCaT.*x(4)^3*x(10)*(x(1)-VCa))/Volm + NPQ*Volmud/Volm*Bmud*x(3)*x(9)*(x(14)-x(15)) + NL*Volmud/Volm*Bmud*x(2)^2*x(8)*(x(13)-x(15)) - Volc/Volm*kPMCA*x(15) + Volc/Volm*Bm*(x(15)-x(16))));
(f*(Bm*(x(15)-x(16)) + pleak*(x(17)-x(16)) - kSERCA*x(16)));
(-f*VolcDIVVoler*(pleak*(x(17)-x(16)) - kSERCA*x(16)))
];
% ODE with the initial condition
xinit = [-65 zeros(1,16)];
t_start = tic;
tspan_ode = [0 time(end)];
[t1, sol] = ode45(f1, time, xinit);
t_end = toc(t_start);
fprintf('\n     [ OK ]');
fprintf('\n     Elapsed %.4f s',t_end);
%% Control Signal & Plasma Membrane Potential
i_stimulation = controlSignal(t1,amp,time,step,npulse,offset);
Vm = sol(:, 1);
figure,
subplot(2,2,[1,2]), plot(t1,i_stimulation, 'LineWidth', 1); hold on; grid minor;
title('Control signal', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$i_m$ [$\mu$A/cm$^2$]', 'Interpreter', 'latex', ...
    'FontSize', 14); 
subplot(2,2,[3,4]), plot(t1, Vm, 'LineWidth', 1); hold on; grid minor;
title('Plasma membrane potential', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$v_m$ [mV]', 'Interpreter', 'latex', ...
    'FontSize', 14); 
%% Activation and Inactivation Functions from ODEs
% Activation Functions
mCaL = sol(:,2); 
mCaPQ = sol(:,3);
mCaT = sol(:,4);
mNa = sol(:,5);
mK = sol(:,6);
mKa = sol(:,7); 
figure, 
subplot(4,1,1), plot(t1, mCaL, 'LineWidth', 1); grid minor;
title('Activation Functions', 'Interpreter', 'latex', 'FontSize', 14);
legm1 = legend('$m_{\mathrm{Ca}_L}$');
set(legm1, 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4,1,2), plot(t1, mCaT, 'LineWidth', 1); grid minor;
legm3 = legend('$m_{\mathrm{Ca}_T}$');
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
set(legm3, 'Interpreter', 'latex', 'FontSize', 14);
subplot(4,1,3), plot(t1, mNa, 'LineWidth', 1); grid minor;
legm4 = legend('$m_{\mathrm{Na}}$');
set(legm4, 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
subplot(4,1,4), plot(t1, mK, 'LineWidth', 1); grid minor;
legm5 = legend('$m_{\mathrm{K}}$');
set(legm5, 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
% % Inactivation Functions
hCaL = sol(:,8);
hCaPQ = sol(:,9);
hCaT = sol(:,10);
hNa = sol(:,11);
hKa = sol(:,12);
figure, 
subplot(3,1,1), plot(t1, hCaL, 'LineWidth', 1); grid minor;
title('Inactivation Functions', 'Interpreter', 'latex', 'FontSize', 14);
legh1 = legend('$h_{\mathrm{Ca}_L}$');
set(legh1, 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,1,2), plot(t1, hCaT, 'LineWidth', 1); grid minor;
legh3 = legend('$h_{\mathrm{Ca}_T}$');
set(legh3, 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,1,3), plot(t1, hNa, 'LineWidth', 1); grid minor;
legh4 = legend('$h_{\mathrm{Na}}$');
set(legh4, 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
%% Calcium Concentrations from ODEs
CaLo = sol(:,13);
CaPQo = sol(:,14);
Cam = sol(:,15);
Cac = sol(:,16);
Caer = sol(:,17);
figure, 
subplot(2,2,1), plot(t1, CaLo, 'LineWidth', 1); hold on; grid minor;
title('$[\mathrm{Ca}]_{L}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14);
subplot(2,2,2), plot(t1, Cam, 'LineWidth', 1); hold on; grid minor;
title('$[\mathrm{Ca}]_m$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14);
subplot(2,2,3), plot(t1, Cac, 'LineWidth', 1); hold on; grid minor;
title('$[\mathrm{Ca}]_c$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14); 
subplot(2,2,4), plot(t1, Caer, 'LineWidth', 1); hold on; grid minor;
title('$[\mathrm{Ca}]_{ER}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14);
%% Exocytosis Equations
% Release Rates
rateCaL = mCaL.^2.*hCaL.*(CaLo.^nL)./(CaLo.^nL + KL.^nL) + ...
    (1 - mCaL.^2.*hCaL).*(Cam.^nL)./(Cam.^nL + KL.^nL);
rateCam = (Cam.^Km)./(Cam.^nm + Km.^nm);
figure,
subplot(3,1,1), 
plot(t1, rateCaL, 'LineWidth', 1, 'LineStyle', '--'); hold on; grid minor;
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L/ms', 'Interpreter', 'latex', 'FontSize', 14);
title('$\mathcal{R}_{\mathrm{Ca}_L}$', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,1,2), 
plot(t1, rateCam, 'LineWidth', 1, 'LineStyle', '-.'); hold on; grid minor;
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L/ms', 'Interpreter', 'latex', 'FontSize', 14);
title('$\mathcal{R}_{\mathrm{Ca}_m}$', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,1,3), 
plot(t1, (rateCaL+rateCam), 'LineWidth', 2); hold on; grid minor;
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L/ms', 'Interpreter', 'latex', 'FontSize', 14);
title('$\mathcal{R}^{\mathrm{(neuron)}}$', 'Interpreter', 'latex', 'FontSize', 14);
%% Released Concentrations
cExoL = zeros(size(rateCaL));      
cExom = zeros(size(rateCam));
for i = 1:length(rateCaL)
    cExoL(i) = trapz(rateCaL(1:i))*step;
    cExom(i) = trapz(rateCam(1:i))*step;
end
figure,
subplot(3,1,1),
plot(t1, cExoL, 'LineWidth', 1, 'LineStyle', '--'); hold on; grid minor;
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L = $\mu$M', 'Interpreter', 'latex', 'FontSize', 14);
title('$c_{\mathrm{Ca}_L}$', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,1,2), 
plot(t1, cExom, 'LineWidth', 1, 'LineStyle', '-.'); hold on; grid minor;
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L = $\mu$M', 'Interpreter', 'latex', 'FontSize', 14);
title('$c_{\mathrm{Ca}_m}$', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,1,3), 
plot(t1, (cExoL+cExom), 'LineWidth', 2); hold on; grid minor;
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L = $\mu$M', 'Interpreter', 'latex', 'FontSize', 14);
title('$c_{\mathrm{Tx}}^{\mathrm{(neuron)}}$', 'Interpreter', 'latex', 'FontSize', 14);