close all; clear all; clc;

% time setting
T = 0.05;
t = 0:T:20;
N = length(t);

% obstacles
o = [[-1;1], [1;-1]];
D   = 1.0;
D_o = sqrt(2)/2;
D_s = sqrt(2)/1.01;


% initial state
p0     = [-2; -2];
x0     = zeros(8,1);
X  = [p0; x0];

%
H = graphic([], 0);

for controller_idx = 1:3
    X  = [p0; x0];
    for i = 1:N
        [dX, ctrl] = closed_loop(X, t, controller_idx);
        H = graphic(H, t(i));
        X = X + dX*T;
    end
end



function [dX, ctrl] = closed_loop(X, t, controller_idx)
o   = evalin("base", 'o');
D   = evalin("base", 'D');
D_o = evalin("base", 'D_o');
D_s = evalin("base", 'D_s');

A=[-1.5828    2.9188         0         0         0         0         0         0
    -2.9188   -1.5828         0         0         0         0         0         0
    0         0   -2.6833    7.1816         0         0         0         0
    0         0   -7.1816   -2.6833         0         0         0         0
    0         0         0         0   -2.5615    6.8558         0         0
    0         0         0         0   -6.8558   -2.5615         0         0
    0         0         0         0         0         0   -2.1391    3.7051
    0         0         0         0         0         0   -3.7051   -2.1391];

B=[ 1.6527         0
    0.6473         0
    1.4972         0
    0.9178         0
    0    1.5791
    0    0.8422
    0    1.5056
    0   -2.2906];
C=[ 0.7761   -1.9815         0         0    2.1315    2.4144         0         0
    0         0   -2.2000   -2.8153         0         0   -1.5060   -0.9899];

p     = X(1:2,   1);
x     = X(3:10,  1);
v_c   = [5.0; 5.0];

if(controller_idx == 1)
    [v_ref, ~, ~] = QP_with_original_feasible_set(p, o, D, D_o, D_s, v_c);
elseif(controller_idx == 2)
    [v_ref, ~, ~] = QP_with_original_feasible_set(p, o, D, D_o, D_s, v_c);
else
    [v_ref, ~, ~] = QP_with_reshaped_feasible_set(p, o, D, D_o, D_s, v_c);
end

ctrl.v_ref = v_ref;
ctrl.v_c   = v_c;
ctrl.u     = v_ref*0;
ctrl.v     = C*x;

% update
dx      =  A*x + B*v_ref;
v       =  C*x;
dp      =  v;
if(controller_idx == 1)
    dX      = [v_ref; dx];
else
    dX      = [dp; dx];
end

end

function [v_ref, A, b] = QP_with_original_feasible_set(p, o, D, D_o, D_s, vc)

opt    = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'active-set');
n_o    = size(o, 2);
A    = zeros(n_o, 2);
b    = zeros(n_o, 1);

for i_o = 1:n_o
    A(i_o, :) = -(p-o(:, i_o))'/norm(p-o(:, i_o));
    b(i_o, 1) = -alpha_h((1/norm(p-o(:, i_o)) - 1/D_s));
end

[v_ref, ~, exitflag, ~, ~] = quadprog(eye(2), -vc, A, b, [], [], [], [], zeros(2,1), opt);

if(exitflag == -2)
    v_ref = zeros(2,1);
end
end

function [v_ref, A, b] = QP_with_reshaped_feasible_set(p, o, D, D_o, D_s, vc)
persistent A_L;
if(isempty(A_L))
    n_L  = 11;
    A_L  = get_positive_basis(n_L);
end

opt    = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'active-set');
n_o    = size(o, 2);
A    = zeros(n_o, 2);
b    = zeros(n_o, 1);

for i_o = 1:n_o
    A(i_o, :) = -(p-o(:, i_o))'/norm(p-o(:, i_o));
    b(i_o, 1) = -alpha_h((1/norm(p-o(:, i_o)) - 1/D_s));
end
k_phi = 5;
[b_L, u_L]= reshape_b_L(A_L, A, b, k_phi);

[v_ref, ~, exitflag, ~, ~] = quadprog(eye(2), -vc, A_L, b_L, [], [], [], [], zeros(2,1), opt);

if(exitflag == -2)
    v_ref = zeros(2,1);
end
end

function y = alpha_h(s)
k = [1, 2];
if(s >= 0)
    y = k(1)*s;
else
    y = k(2)*s;
end
end

function [b_L, nu_s]= reshape_b_L(A_L, A, b, k_psi)
n_l  = size(A_L, 1);
c_A  = cos(2*pi/n_l);
nu_s = A' * min(b, 0);
psi  = (b-A*nu_s)' .* max(A_L * A', c_A) + k_psi*max(c_A - A_L*A', 0);
b_L  = A_L*nu_s + min(psi, [], 2);
end

function A = get_positive_basis(N)
    theta = linspace(0, 2*pi, N+1)';
    A     = [cos(theta(1:end-1)), sin(theta(1:end-1))];
end



function H = graphic(H, t_now)
persistent ctrl_last_proposed ctrl_last_standard;

o   = evalin("base", 'o');
D   = evalin("base", 'D');
D_o = evalin("base", 'D_o');
D_s = evalin("base", 'D_s');

X   = evalin("base", 'X');
T   = evalin("base", 'T');

p   = X(1:2,   1);
x   = X(3:10,  1);

delta_h = D_s(1)-D(1);

try
    ctrl = evalin("base", 'ctrl');
    controller_idx = evalin("base", 'controller_idx');
catch

end

if(isempty(H))
    clr   = {'k',  'r', 'b'};
    style = {'-.', ':', '-'};

    H.fig1 = figure(1);
    hold on; grid on; axis equal; box on;
    H.axes1 = gca;
    set(gcf, 'position', [50, 50, 650, 600], 'color', 'w');
    set(gca, 'fontsize', 16, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
    axis([-3, 3, -3, 3]);
    
    xlabel('$[p]_1$', 'Interpreter', 'latex');
    ylabel('$[p]_2$', 'Interpreter', 'latex');
    theta = linspace(0,2*pi,36);
    for i = 1:size(o, 2)
        H.f1_obs_margin{i} = patch(o(1,i)+cos(theta)*(D(1)+delta_h), o(2,i)+sin(theta)*(D(1)+delta_h), 'c', 'handlevisibility', 'off', 'displayname', 'Obstacles', 'linestyle', 'none', 'facealpha', 0.2);
        temp = patch(o(1,i)+cos(theta)*D(1), o(2,i)+sin(theta)*D(1), 'k', 'handlevisibility', 'off', 'displayname', 'Obstacles', 'linestyle', '-', 'facealpha', 1.0);
        H.f1_obs{i} = hatchfill(temp, 'single', 45, 10, [1.0, 1.0, 1.0]);
        H.f1_obs{i}.HandleVisibility = 'off';
    end

    H.f1_traj{1} = animatedline('Color',clr{1},'DisplayName','Case 1','LineStyle',style{1}, 'LineWidth', 2);
    H.f1_traj{2} = animatedline('Color',clr{2},'DisplayName','Case 2','LineStyle',style{2}, 'LineWidth', 2);
    H.f1_traj{3} = animatedline('Color',clr{3},'DisplayName','Case 3','LineStyle',style{3}, 'LineWidth', 2);
    legend('Location','northwest');

    H.fig5  = figure(5);
    set(gcf, 'position', [600, 100, 900, 500], 'color', 'w');
    subplot(2,1,1);
    H.axes5_1 = gca;
    hold on; grid on; box on;
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$[v^*]_1$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
    H.f5_vref1_1 = animatedline('Color',clr{1},'DisplayName','Case 1','LineStyle',style{1}, 'LineWidth', 2);
    H.f5_vref1_2 = animatedline('Color',clr{2},'DisplayName','Case 2','LineStyle',style{2}, 'LineWidth', 2);
    H.f5_vref1_3 = animatedline('Color',clr{3},'DisplayName','Case 3','LineStyle',style{3}, 'LineWidth', 2);
    legend('Location','northwest');

    subplot(2,1,2);
    H.axes5_2 = gca;
    hold on; grid on; box on;
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$[v^*]_2$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
    H.f5_vref2_1 = animatedline('Color',clr{1},'DisplayName','Case 1','LineStyle',style{1}, 'LineWidth', 2);
    H.f5_vref2_2 = animatedline('Color',clr{2},'DisplayName','Case 2','LineStyle',style{2}, 'LineWidth', 2);
    H.f5_vref2_3 = animatedline('Color',clr{3},'DisplayName','Case 3','LineStyle',style{3}, 'LineWidth', 2);
    legend('Location','northwest');


else

    addpoints(H.f1_traj{controller_idx}, p(1), p(2));

    if(controller_idx == 1)
        addpoints(H.f5_vref1_1, t_now, ctrl.v_ref(1));
        addpoints(H.f5_vref2_1, t_now, ctrl.v_ref(2));
    elseif(controller_idx == 2)
        addpoints(H.f5_vref1_2, t_now, ctrl.v_ref(1));
        addpoints(H.f5_vref2_2, t_now, ctrl.v_ref(2));
    elseif(controller_idx == 3)
        addpoints(H.f5_vref1_3, t_now, ctrl.v_ref(1));
        addpoints(H.f5_vref2_3, t_now, ctrl.v_ref(2));
    end
end

drawnow limitrate;
end


