% YAO ZHANG 16-11-2022

clear all
close all
clc

%% Single point absorber WEC parameters:
Ts = 0.1;
g = 9.81;
r = 0.35;  % Raidus;
rho = 1025; %  The average density of seawater at the ocean surface;
m =242+83.5;
ks = pi*r^2*rho*9.8; % Stifness;
D = 2e5;
Dv = 0;

Df = 0
Ptime = 50;
Prediction = Ptime/Ts;

% 10nd-order modeling of WEC 
% radiation force
Ar= [];
Br = [];
Cr = [];
% [Ar Br Cr Dr] = tf2ss(num,dem)
%excitation force
[n_r,~] = size(Ar)
Ae= [ ];
Be = [];
Ce = [];
[n_e,~] = size(Ae)

%% WEC state-space continous time

Ac = [0 1 zeros(1,n_r+n_e)
    -ks/m -Df/m Cr/m -Ce/m
    zeros(n_r,1) Br Ar zeros(n_r,n_e)
    zeros(n_e,2+n_r) Ae]
Buc = [0 ;1000/m;zeros(n_r+n_e,1)]
Bwc = [0;0;zeros(n_r,1);Be]
Cc = [1 0 zeros(1,n_e+n_r)]
Cz = [0 1 zeros(1,n_e+n_r)]
[nx,~] =size(Ac)
%% WEC state-space discrete time
H = ss(Ac, [Bwc, Buc], eye(nx), zeros(nx,2));
Hd = c2d(H, Ts,'zoh');
[A, B, ~, ~] = ssdata(Hd);
Bw = B(:,1);
Bu = B(:,2);

%% wave data
N=50;
tend = 200;



% Interpolation of wave amplitud ;
ti = (0:Ts:(trunc-1)*TS)'; 
w=interp1(time,w_origin,ti);
dim_wi = size(w,1);

% Wave speed:
wd_origin = (w_origin(2:trunc)-w_origin(1:trunc-1))/TS;
ts = (0:Ts:(trunc-1)*TS-N*Ts)';
times = wdata(1:trunc-1,1);
wd = interp1(times,wd_origin,ts);

Nn=tend/Ts+1;
Nm = Nn;
figure(1)
plot(ti(1:Nn), w(1:Nn),'m-','LineWidth',1.5)
xlabel('Time (s)','FontSize',14,'FontName','Times New Roman')
ylabel('Amplitude (m)','FontSize',14,'FontName','Times New Roman')
title('wave prediction')
grid on;

%% Controller design 
plotstyle = 'c--';
control_horizon = 20;
% controller parameter tuning 
R = ;
T = ;
q1 = ;
q2 = ;

Q = blkdiag(q1,q2,eps*eye(nx-2));

S = dare(A,Bu,Q,R,T*Cz',eye(nx));
Ks=-(Bu'*S*Bu+R)\Bu';
Kw=-Ks*S*Bw;
Kx=-(Bu'*S*Bu+R)\(Bu'*S*A+T*Cz);
Phi = (A+Bu*Kx)';   
Psi = zeros(nx,control_horizon);
for i = 1:control_horizon
Psi(:,i) = Phi^(i-1)*S*Bw;
end
Kd = Ks*Psi;


x_loc = zeros(10,1);
PE_loc = 0;
x_cau = zeros(10,1);
PE_cau = 0;
t0 = clock;

for i=1:Nm
% Noncausal linear control
u_loc = Kx*x_loc+Kd*w(i:i+control_horizon-1);
U_loc(i) = u_loc;

y_loc= x_loc(2);
Y_loc(i) = y_loc;

P_loc(i)=-u_loc*y_loc*Ts;
PE_loc=PE_loc+P_loc(i);
E_loc(i)=PE_loc;

x_loc=A*x_loc+Bu*u_loc+Bw*w(i);

% Noncausal control
u_cau = Kx*x_cau;
U_cau(i) = u_cau;

y_cau= x_cau(2);
Y_cau(i) = y_cau;

P_cau(i)=-u_cau*y_cau*Ts;
PE_cau=PE_cau+P_cau(i);
E_cau(i)=PE_cau;

x_cau=A*x_cau+Bu*u_cau+Bw*w(i);
end

disp('The noncausal control energy generated is:')
disp(max(E_loc))



%% simulation
figure(2)
plot(ts(1:Nm),E_loc,'b',ts(1:Nm),E_cau,'k','LineWidth',1.5)
legend('Noncausal','Causal')
title('Noncausal vs Causal')
xlabel('Time (s)')
ylabel('Energy (kJ)')
figure(3)
plot(ts(1:Nm),U_loc,'b',ts(1:Nm),U_cau,'k','LineWidth',1.5)
legend('Noncausal','Causal')
title('Noncausal vs Causal')
xlabel('Time (s)')
ylabel('PTO force (N)')
