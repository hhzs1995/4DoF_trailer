%% Nonlinear simulations for two-wheeled trailers
% 4 DoF model, q = [psi, theta, fi, u]

% based on: H. Z. Horvath, D. Takacs, M. Antali: 
% Effects of differential braking and tire contact modeling 
% on the stability and nonlinear dynamics of towed two-wheeled trailers 
% (2026). Manuscript submitted for publication.

clear all;
close all;
clc

% options
t_end = 20; % end of simulation
y0 = [0 0 0 0 -0.5 0 0 0]; % psi'(0) as initial condition
h_max = 0.001; % max. step size
tau = 0.01; % time delay, s
options_dde = ddeset('InitialY', y0,'MaxStep',h_max,'RelTol',1e-5,'AbsTol',1e-7);

% simulation with dde23
sol = dde23(@trailer_4DoF_goveq, tau, @trailer_hist,[0, t_end], options_dde);

% plotting
figure
% psi(t)
subplot(4,1,1)
plot(sol.x,sol.y(1,:),'k','Linewidth', 1)
hold on;
line([0 t_end],[0 0],'Color','black','LineStyle','--')
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 10;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$\psi [\mathrm{rad}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
% theta(t)
subplot(4,1,2)
plot(sol.x,sol.y(2,:),'k','Linewidth', 1)
hold on;
line([0 t_end],[0 0],'Color','black','LineStyle','--')
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 10;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$\vartheta [\mathrm{rad}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
% phi(t)
subplot(4,1,3)
plot(sol.x,sol.y(3,:),'k','Linewidth', 1)
hold on;
line([0 t_end],[0 0],'Color','black','LineStyle','--')
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 10;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$\varphi [\mathrm{rad}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
% u(t)
subplot(4,1,4)
plot(sol.x,sol.y(4,:),'k','Linewidth', 1)
hold on;
line([0 t_end],[0 0],'Color','black','LineStyle','--')
ax = gca;
ax.FontName = 'TimesNewRoman';
ax.FontSize = 10;
xlabel('$t [\mathrm{s}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');
ylabel('$u [\mathrm{m}]$','Interpreter','latex','FontSize',8,'FontName','TimesNewRoman');

function s = trailer_hist(t)

    s = zeros(8,1); 
    
end