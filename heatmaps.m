% PARAMS
dt = 0.1;
g = -10.;
l = 1.; % length pendulum.
m = 1.; 
beta = 1.; % drag coeff. 
A = [1. dt; -dt*g/l 1-beta/m*dt]; % evolution matrice, x(k+1) = A x(k) + U. 
k1 = 1; % num of forward steps for the supervisor.
k2 = 25; % num of backward steps. 
max_u = 6; % max input. 
eta = 0; % noise on the policy: a percentage of the max interval. 
% END PARAMS

close all;

figure;

n = 80; % dimension of the heatmap. 

h1 = stochastic_heatmap(A, dt, max_u, n, eta, 25);
h2 = supervision_heatmap(A, dt, max_u, k1, k2, n, eta, 25);
% save('super200.mat','h2');
h = heatmap(h1, 'Colormap', winter, 'GridVisible','off');
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
colorbar('off')

% Plot contour line related to supervised set. 

x = 1:n;
y = 1:n;
[X,Y] = meshgrid(x,y);
Z = ones(n);
for i = 1:n
    for j = 1:n   
        Z(i, j) = h2(X(i, j), Y(i, j));
    end
end
ax1 = axes;
Z = rot90(Z);
[x1,y1] = contour(ax1, X, Y, Z, [1. 1.],'red', 'linewidth', 2.);
ax1.Color = 'none';
set(gca,'Ytick',[]) 
set(gca,'Xtick',[])

% Plot backward set.

ax2 = axes;
hold(ax2,'on')
T = zonotope(interval([0.; 0.],[0.; 0.]));
U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
bX = k_step_backward(T, U, A, dt, k2);
plot(bX, [1,2], 'k', 'linewidth', 2.);
ax2.Color = 'none';
xlim([-1 1]);
ylim([-2 2]);

% Plot test points. 

% ax3 = axes;
% hold(ax2,'on')
% scatter(starting_pos(1, :), starting_pos(2, :), 20, 'm', 'filled');
% ax3.Color = 'none';
% xlim([-1 1]);
% ylim([-2 2]);

function h = supervision_heatmap(A, dt, max_u, k1, k2, dim, eta, horizon)
    heat = ones(dim);
    X1 = linspace(-1, 1, dim);
    X2 = linspace(-2, 2, dim);
    T = zonotope(interval([0.; 0.],[0.; 0.]));
    U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
    B = k_step_backward(T, U, A, dt, k2);
    for i = 1:dim
        for j = 1:dim
            x2 = X2(dim + 1 - i);
            x1 = X1(j);
            x = [x1; x2];
            for t = 1:horizon
                if mod(t, k1) == 0
                    %X = zonotope(interval([x(1,1); x(2, 1)],[x(1,1); x(2, 1)]));
                    W = get_max_w_quick(x, B, U, A, dt);
                 end
                 u = pd_control(x, max_u, eta);
                 if not(isempty(W.vertices))
                     u = supervision(u, W);
                 end
                x = move(x, u, A, dt);
                if abs(x(1, 1)) > 1.
                    break;
                end
            end
                heat(i, j) = norm(x);
        end
    end
    h = heat;
end

function h = stochastic_heatmap(A, dt, max_u, dim, eta, horizon)
    heat = ones(dim);
    X1 = linspace(-1, 1, dim);
    X2 = linspace(-2, 2, dim);
    for i = 1:dim
        for j = 1:dim
            x2 = X2(dim + 1 - i);
            x1 = X1(j);
            x = [x1; x2];
            for t = 1:horizon
                u = pd_control(x, max_u, eta);
                %u = pd_control(x, max_u);
                x = move(x, u, A, dt);
                if abs(x(1, 1)) > 1.
                    break;
                end
            end
            heat(i, j) = norm(x);
        end
    end
    h = heat;
end

function res = supervision(u, w)
    inter = w.vertices;
    sup = inter(2, 2);
    low = inter(2, 1);
    if u < low
        u = low;
    elseif u > sup
        u = sup;
    end
    res = u;
end

function res = k_step_forward(x, u, A, dt, k)
    for i = 1:k
        x = A*x + dt*u;
    end
    res = x;
end
function res = k_step_backward(x, u, A, dt, k)
    for i = 1:k
        x = inv(A) *  (x + -dt*u);
    end
    res = x;
end

function [w, w_box] = get_max_w_quick(X, B, U, A, dt)
    w_box = polytope((1/dt)*(-A*X + B));
    w = w_box & U;
end

function res = clip(x, inf , sup)
    if x < inf
        x = inf;
    elseif x > sup
        x = sup;
    end
    res = x;
end


function res = pd_control(x, max_u, eta)
    g = -10.;
    l = 1.;
    m = 1.;
    beta = 1.;
    a = -g/l;
    b = -beta/m;
    u = -(a+1)*x(1,1) - (b+2)*x(2, 1);
    res = clip(u + eta * (2 * max_u *rand() - max_u), -max_u, max_u);
end

function res = move(x, u, A, dt)
        res = A*x + [0; dt]*u;
end

% function res = move(x, u, A, dt)
%     A = [0 1; 10 -1]; 
%     inva = inv(A);
%     K = x + inva * [0; u];
%     x  = expm(dt * A) * K - inva * [0; u]; 
%     res = x; 
% end