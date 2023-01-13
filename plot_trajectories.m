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
horizon = 25; 

% END PARAMS


n_points = 60;
cc = hsv(n_points);
points = [2 * rand(1, n_points) - 1; 4 * rand(1, n_points) - 2.]; % initial conditions.
T = [0; 0]; 
U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
B = k_step_backward(T, U, A, dt, k2);
close all;

figure;
hold on;
xlim([-1 1]);
ylim([-2 2]);

plot(bX, [1,2], 'k', 'linewidth', 2.); % Plot backward set

for i = 1:n_points
    x = points(:, i);
    traj = ones(2, horizon+1);
    traj(:, 1) = x;
    col = [.9 .9 .9];
    for t = 1:horizon
        if mod(t, k1) == 0
            W = get_max_w_quick(x, B, U, A, dt);
         end
         u = pd_control(x, max_u, eta);
         if not(W.inf == -Inf)
             [u, super] = supervision(u, W);
         end
         if super
            col = [0. .5 0.];
            scatter(x(1, 1), x(2, 1), 80, 'r', '*'); 
         end
         
        x = move(x, u, A, dt);
        traj(:, t+1) = x;  
    end
    plot(traj(1, :), traj(2, :),'color', col, 'linewidth', 1.);
    scatter(traj(1, 1), traj(2, 1),[], col, '>', 'filled'); 
end 
 
% Plot contour 
x = 1:n;
y = 1:n;
[X,Y] = meshgrid(x,y);
Z = ones(n);
for i = 1:n
    for j = 1:n   
        Z(i, j) = h1(X(i, j), Y(i, j));
    end
end
ax1 = axes;
Z = rot90(Z);
[x1,y1] = contour(ax1, X, Y, Z, [1. 1.],'b--', 'linewidth', 2.);
ax1.Color = 'none';
set(gca,'Ytick',[]) 
set(gca,'Xtick',[])


function [res, super] = supervision(u, w)
    sup = w.sup;
    low = w.inf;
    super = 0;
    if u < low
        u = low;
        super = 1;
    elseif u > sup
        u = sup;
        super = 1; 
    end
    res = u;
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
    temp = w.interval;
    w = temp(2, 1);
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

function res = move(x, u, A, dt)
        res = A*x + [0; dt]*u;
end
