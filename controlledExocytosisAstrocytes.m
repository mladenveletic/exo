%% For Abel cluster
% pc = parcluster
% pc.JobStorageLocation = getenv('SCRATCH')
% parpool(pc,12)
%% For local computation
% delete(gcp('nocreate'))
% parpool('local',4)
%%
step = 1;
npulse = 1;
time = 0:step:50e3;  %[ms]
amp = 10;            %[mV]
IP0 = 0.160;         %[muM]
offset = 12500;
%% Parameter Set: ---
% ASTROCYTE CELLS
% Nominal values for electrical activity
% Voltages given in [mV]
% Conductances given in [nS or, equivalently, mS/cm^3]
% Time constants given in [ms]
% Capacitance given in [\muF/cm^2]
% Parameters from alpha cells
VCa = 65;
gCaL = 0.7;
gCaPQ = 0.6;
gCaT = 0.4;
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
VtaumK = -10;
StaumK = 25;
% Nominal values for calcium dynamics and exocytosis
f = 0.01;
Volmud = 2.618e-19; % [L]
Volc = 5.725e-13;   % [L]
NPQ = 200;
Bmud = 264;         % [1/ms]
kPMCA = 0.3;        % [1/ms]
nPQ = 4;
nL = 4;
nm = 4;
alpha = 5.18e-15;   % [\mumol/pA/ms]
Volm = 5.149e-14;   % [L]
NL = 200;
Bm = 0.128;         % [1/ms]
KPQ = 2;
KL = 50;            % [muM]
Km = 2;             % [muM]
% Nominal values for calcium dynamics in astrocytes
MS = 1e3;
c0 = 2;             % [muM]
c1 = 0.185;
v1 = 6/MS;          % [1/ms]
v2 = 0.11/MS;       % [1/ms]
v3 = 0.9/MS;        % [muM/ms]
k3 = 0.1;           % [muM]
d1 = 0.13;
d2 = 1.049;
d3 = 0.943;
d5 = 0.082;
a2 = 0.5/MS;        % [1/(muM ms)]
r_IP = 0.04/MS;        % [muM/ms]
tau_IP = 1/0.000140;% [ms]
%% Derivatives
% x(1) = mCaL;
% x(2) = mCaPQ;
% x(3) = mCaT;
% x(4) = hCaL;
% x(5) = hCaPQ;
% x(6) = hCaT;
% x(7) = CaL;
% x(8) = CaPQ;
% x(9) = Cam;
% x(10) = Cac;  -->modified according to Li-Rinzel model
% x(11) = Caer; -->modified according to Li-Rinzel model
% x(12) = IP3;  -->added according to Li-Rinzel model
% x(13) = hIP3; -->added according to Li-Rinzel model
%--------------------------------------------------------------------------
pot = -70; % Resting potential
f1 = @(t,x) [
    (funcminf(pot+controlSignal(t,amp,time,step,npulse,offset), VmCaL, SmCaL) - x(1))./functaum(pot+controlSignal(t,amp,time,step,npulse,offset), taumVCaL, VtaumCaL, StaumCaL, taum0CaL);
    (funcminf(pot+controlSignal(t,amp,time,step,npulse,offset), VmCaPQ, SmCaPQ) - x(2))./functaum(pot+controlSignal(t,amp,time,step,npulse,offset), taumVCaPQ, VtaumCaPQ, StaumCaPQ, taum0CaPQ);
    (funcminf(pot+controlSignal(t,amp,time,step,npulse,offset), VmCaT, SmCaT) - x(3))./functaum(pot+controlSignal(t,amp,time,step,npulse,offset), taumVCaT, VtaumCaT, StaumCaT, taum0CaT);
    (funcminf(pot+controlSignal(t,amp,time,step,npulse,offset), VhCaL, ShCaL) - x(4))./functaum(pot+controlSignal(t,amp,time,step,npulse,offset), tauhVCaL, VtauhCaL, StauhCaL, tauh0CaL);
    (funcminf(pot+controlSignal(t,amp,time,step,npulse,offset), VhCaPQ, ShCaPQ) - x(5))./functaum(pot+controlSignal(t,amp,time,step,npulse,offset), tauhVCaPQ, VtauhCaPQ, StauhCaPQ, tauh0CaPQ);
    (funcminf(pot+controlSignal(t,amp,time,step,npulse,offset), VhCaT, ShCaT) - x(6))./functaum(pot+controlSignal(t,amp,time,step,npulse,offset), tauhVCaT, VtauhCaT, StauhCaT, tauh0CaT);
    (-f*alpha*gCaL.*(pot+controlSignal(t,amp,time,step,npulse,offset)-VCa)./NL./Volmud - f*Bmud*(x(7)-x(9)));
    (-f*alpha*gCaPQ.*(pot+controlSignal(t,amp,time,step,npulse,offset)-VCa)./NPQ./Volmud - f*Bmud*(x(8)-x(9)));
    (f*(-alpha*(gCaT.*x(3)^3*x(6)*(pot+controlSignal(t,amp,time,step,npulse,offset)-VCa))/Volm + NPQ*Volmud/Volm*Bmud*x(2)*x(5)*(x(8)-x(9)) + NL*Volmud/Volm*Bmud*x(1)^2*x(4)*(x(7)-x(9)) - Volc/Volm*kPMCA*x(9) + Volc/Volm*Bm*(x(9)-x(10))));
    -c1*v1*(x(12)/(x(12)+d1)*x(10)/(x(10)+d5))^3*x(13)^3*(x(10)-x(11)) - c1*v2*(x(10)-x(11)) - v3*x(10)^2/(k3^2+x(10)^2);
    -1/c1*(-c1*v1*(x(12)/(x(12)+d1)*x(10)/(x(10)+d5))^3*x(13)^3*(x(10)-x(11)) - c1*v2*(x(10)-x(11)) - v3*x(10)^2/(k3^2+x(10)^2));
    (IP0 - x(12))/tau_IP + r_IP*controlSignal(t,amp,time,step,npulse,offset);
    a2*d2*(x(12)+d1)/(x(12)+d3)*(1-x(13)) - a2*x(10)*x(13);
    ];
% ODE with initial condition
xinit = [zeros(1,10) c0/c1 IP0 0];
t_start = tic;
% Tim Franklin (2020). ODE Progress Bar and Interrupt (https://www.mathworks.com/matlabcentral/fileexchange/9904-ode-progress-bar-and-interrupt), MATLAB Central File Exchange. Retrieved January 16, 2020.
options=odeset('RelTol',1e-3,'AbsTol',1e-4,'InitialStep',0.001,'MaxStep',30,...
    'OutputFcn',@odeprog,'Events',@odeabort);
[t1, sol] = ode15s(f1, time, xinit, options);
t_end = toc(t_start);
fprintf('\n     [ OK ]');
fprintf('\n     Elapsed %.4f s',t_end)
%% Control Signal & Plasma Membrane Potential
Vm = controlSignal(t1,amp,time,step,npulse,offset);
figure,
subplot(2,1,1), plot(t1, Vm, 'LineWidth', 1); hold on, grid minor,
title('Control signal', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$v_m$ [mV]', 'Interpreter', 'latex', 'FontSize', 14);
IP3 = sol(:,12);
subplot(2,1,2), plot(t1, IP3, 'LineWidth', 1); hold on, grid minor,
title('$[\mathrm{IP}_3]$ concentration', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$[\mathrm{IP}_3]$ [$\mu$mol/L]', 'Interpreter', 'latex', 'FontSize', 14);
%% Activation and Inactivation Functions from ODEs
% Activation Functions
mCaL = sol(:,1);
mCaPQ = sol(:,2);
mCaT = sol(:,3);
mIP3inf = sol(:,12)./(sol(:,12)+d1).*sol(:,10)./(sol(:,10)+d5);
% Inactivation Functions
hCaL = sol(:,4);
hCaPQ = sol(:,5);
hCaT = sol(:,6);
hIP3 = sol(:,13);
%% Calcium Concentrations from ODEs
CaLo = sol(:,7);
CaPQo = sol(:,8);
Cam = sol(:,9);
Cac = sol(:,10);
Caer = sol(:,11);
figure,
subplot(3,2,1), plot(t1, CaLo, 'LineWidth', 1); hold on, grid minor,
title('$[\mathrm{Ca}]_{L}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,2,2), plot(t1, CaPQo, 'LineWidth', 1); hold on, grid minor,
title('$[\mathrm{Ca}]_{N}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,2,3), plot(t1, Cam, 'LineWidth', 1); hold on, grid minor,
title('$[\mathrm{Ca}]_m$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,2,4), plot(t1, Cac, 'LineWidth', 1); hold on, grid minor,
title('$[\mathrm{Ca}]_c$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14);
subplot(3,2,5), plot(t1, Caer, 'LineWidth', 1); hold on, grid minor,
title('$[\mathrm{Ca}]_{ER}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L', 'Interpreter', 'latex', 'FontSize', 14);
%% Exocytosis Equations
% Release Rates
rateCaL = mCaL.^2.*hCaL.*(CaLo.^nL)./(CaLo.^nL + KL.^nL) + ...
    (1 - mCaL.^2.*hCaL).*(Cam.^nL)./(Cam.^nL + KL.^nL);
rateCaPQ = mCaPQ.*hCaPQ.*(CaPQo.^nPQ)./(CaPQo.^nPQ + KPQ.^nPQ) + ...
    (1 - mCaPQ.*hCaPQ).*(Cam.^nPQ)./(Cam.^nPQ + KPQ.^nPQ);
rateCam = (Cam.^Km)./(Cam.^nm + Km.^nm);
figure,
subplot(2,2,1),
plot(t1, rateCaL, 'LineWidth', 1, 'LineStyle', '--'); hold on, grid minor,
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L/ms', 'Interpreter', 'latex', 'FontSize', 14);
title('$\mathcal{R}_{\mathrm{Ca}_L}$','Interpreter', 'latex','FontSize', 14);
subplot(2,2,2),
plot(t1, rateCaPQ, 'LineWidth', 1, 'LineStyle', '-'); hold on, grid minor,
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L/ms', 'Interpreter', 'latex', 'FontSize', 14);
title('$\mathcal{R}_{\mathrm{Ca}_{N}}$','Interpreter', 'latex','FontSize', 14);
subplot(2,2,3),
plot(t1, rateCam, 'LineWidth', 1, 'LineStyle', '-.'); hold on, grid minor,
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L/ms', 'Interpreter', 'latex', 'FontSize', 14);
title('$\mathcal{R}_{\mathrm{Ca}_m}$','Interpreter', 'latex','FontSize', 14);
subplot(2,2,4),
plot(t1, (rateCaL+rateCaPQ+rateCam), 'LineWidth', 2); hold on, grid minor,
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
title('$\mathcal{R}^{\mathrm{(astro)}}$', 'Interpreter', 'latex','FontSize', 14);
%% Released Concentrations
cExoL = zeros(size(rateCaL));
cExoPQ = zeros(size(rateCaPQ));
cExom = zeros(size(rateCam));
for i = 1:length(rateCaL)
    cExoL(i) = trapz(rateCaL(1:i))*step;
    cExoPQ(i) = trapz(rateCaPQ(1:i))*step;
    cExom(i) = trapz(rateCam(1:i))*step;
end
figure,
subplot(2,2,1),
plot(t1, cExoL, 'LineWidth', 1, 'LineStyle', '--'); hold on, grid minor,
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L = $\mu$M', 'Interpreter', 'latex', 'FontSize', 14);
title('$c_{\mathrm{Ca}_L}$','Interpreter', 'latex','FontSize', 14);
subplot(2,2,2),
plot(t1, cExoPQ, 'LineWidth', 1, 'LineStyle', '-'); hold on, grid minor,
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L = $\mu$M', 'Interpreter', 'latex', 'FontSize', 14);
title('$c_{\mathrm{Ca}_{N}}$','Interpreter', 'latex','FontSize', 14);
subplot(2,2,3),
plot(t1, cExom, 'LineWidth', 1, 'LineStyle', '-.'); hold on, grid minor,
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\mu$mol/L = $\mu$M', 'Interpreter', 'latex', 'FontSize', 14);
title('$c_{\mathrm{Ca}_m}$','Interpreter', 'latex','FontSize', 14);
subplot(2,2,4)
plot(t1, (cExoL+cExoPQ+cExom), 'LineWidth', 2); hold on, grid minor,
xlabel('Time [ms]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Relative Concentration', 'Interpreter', 'latex', 'FontSize', 14);
title('$c_{\mathrm{Tx}}^{\mathrm{(astro)}}$', 'Interpreter', 'latex','FontSize', 14);
