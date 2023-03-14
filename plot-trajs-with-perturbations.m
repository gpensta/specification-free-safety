% PARAMS
dt = 0.1;
g = -10.;
l = 1.; % length pendulum.
m = 1.; 
beta = 1.; % drag coeff. 
A = [1. dt; -dt*g/l 1-beta/m*dt]; % evolution matrice, x(k+1) = A x(k) + U. 
k1 = 1; % num of forward steps for the supervisor.
k2 = 25; % num of backward
max_u = 6; % max input. 
eta = 0.03; % noise on the policy: a percentage of the max interval. 
horizon = 25; 


%   END PARAMS
n_points = 1;
cc = hsv(n_points);
%  points = [-0.2; -1];
points = [2 * rand(1, n_points) - 1; 4 * rand(1, n_points) - 2.]; % initial conditions.

T =  zonotope(interval([-0.1; -0.1;],[0.1 ; 0.1])); 
U = dt * max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
W = dt * eta * max_u * zonotope(interval([-1.; -1.],[1. ; 1.]));
B = k_robust_backward(T, U, W, A, k2);
 
Bw =  k_step_backward(T, U+W, A, k2);
close all; figure; hold on; xlim([-1 1]); ylim([-2 2]);

plot(B, [1,2], 'k', 'linewidth', 2.); % Plot backward set
plot(Bw, [1,2], 'k--', 'linewidth', 2.); % Plot backward set
% plot(Br, [1,2], 'k--', 'linewidth', .5); % Plot backward set
% plot([-1:.1:1], -6 - 11 * [-1:.1:1])

for i = 1:n_points
    x = points(:, i);
    traj = ones(2, horizon+1);
    traj(:, 1) = x;
    col = [.9 .9 .9];
    super = 0;
    had_super = 0;
    for t = 1:horizon
        if mod(t, k1) == 0
            Us = get_max_w_quick(x, B, U, A, dt);
         end
         u = pd_control(x, max_u, eta);
         if not(isempty(Us))
             [u, super] = supervision(u, Us);
         end
         if super
            had_super = 1; 
            col = [0. .5 0.];
            scatter(x(1, 1), x(2, 1), 80, 'r', '*'); 
         end
         
        x = move(x, u, A, dt, eta * max_u);
        traj(:, t+1) = x;  
    end
    cols = [0.9 .9 0.9; 0.5 0. 0.;  0.93 0.69 0.13];
    if abs(traj(1, horizon)) > 0.5
        c = 2;
    else
        c = 1;
    end
    
    if had_super
        
        s = 100;
        c = 3;
    else
        s = 25;
    end
    plot(traj(1, :), traj(2, :),'color', col, 'linewidth', 1.);
    scatter(traj(1, 1), traj(2, 1),s, cols(c, :), 'o', 'filled'); 
end 
 
% Plot contour 
n = 100;
x = 1:n;
y = 1:n;
h1 = stochastic_heatmap(A, dt, max_u, n, eta, horizon);
[X,Y] = meshgrid(x,y);
Z = ones(n);
for i = 1:n
    for j = 1:n   
        Z(i, j) = h1(X(i, j), Y(i, j));
    end
end
% ax1 = axes;
% Z = rot90(Z);
% [x1,y1] = contour(ax1, X, Y, Z, [1. 1.],'b--', 'linewidth', 2.);
% ax1.Color = 'none';
% set(gca,'Ytick',[]) 
% set(gca,'Xtick',[])


function [res, super] = supervision(u, w)
    w = w.vertices;
    sup = w(2, 2);
    low = w(2, 1);
    range = sup - low;
    super = 0;
    if u < low
        u = low; %+ .5 * range;
        super = 1;
    elseif u > sup
        u = sup; % - .5 * range;
        super = 1; 
    end
    res = u;
end
function res = k_step_backward(X, U, A, k)
    for i = 1:k
        X = inv(A) *  (X + -1*U);
        plot(X, [1,2], 'b', 'linewidth', 1.,'FaceColor', 'none');
    end
    res = X;
end

function res = k_robust_backward(X, U, W, A, k)
    for i = 1:k
        d1 = minkDiff(X, W);
        if isempty(d1) 
           break 
        else
            temp = d1 + -1*U;
            d2 = minkDiff(temp, W);
            if isempty(d2)
                break
            else
                X = inv(A) * d2;
                plot(X, [1,2], 'm', 'linewidth', 1.,'FaceColor', 'none');
            end
        end
    end
    res = X;
end


function [w, w_box] = get_max_w_quick(X, B, U, A, dt)
    w_box = mptPolytope((-A*X + B));
%     w_box = mptPolytope(w_box);
    w = (1/dt)  * (w_box & U);
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
    res = clip(u, -max_u, max_u);
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
                x = move(x, u, A, dt, eta * max_u);
                if abs(x(1, 1)) > 1.
                    break;
                end
            end
            heat(i, j) = norm(x);
        end
    end
    h = heat;
end

function res = move(x, u, A, dt, max_w)
        res = A*x + [0; dt]*u + (2 * max_w *rand() - max_w)*dt;
end

function y = map(x, a, b, c, d)
    y = c + ((d - c)*(x-a))/(b -a);
end
