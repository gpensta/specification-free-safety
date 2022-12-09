% PARAMS
dt = 0.1;
g = -10.;
l = 1.;
m = 1.;
beta = 1.;
A = [1. dt; -dt*g/l 1-beta/m*dt];
k1 = 10;
k2 = 9;
max_u = 6;
% END PARAMS

% close all

% figure;
% own = [winter; 1. 1. 1.];
% heatmap(h1c, 'ColorLimits',[0 3], 'Colormap', own); 

horizon = 100;
eta = 0.;

% (26, 48) (15, 33) (26, 1) (19, 24)

dx = 2 / 50;
dy = 4 / 50;

starting_pos = [-1+26*dx -1+15*dx -1+26*dx -1+19*dx -1+30*dx -0.5; 2-48*dy 2-33*dy 2-1*dy 2-24*dy 2-45*dy -0.5]; 
f = figure;
for k =  1:length(starting_pos) 
    x = [starting_pos(1, k); starting_pos(2, k)];
    y = x;

    naive_traj = ones(2, horizon);
    naive_u = ones(1, horizon);

    for t=1:1:horizon
        u = stochastic_control(x, max_u, eta);
        naive_u(t) = u;
        x = move_exact(x, u, dt); 
        naive_traj(1, t) = x(1, 1);
        naive_traj(2, t) = x(2, 1);
    end

    x = y;
    supervised_traj = ones(2, horizon);
    supervised_u = ones(3, horizon);
    U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
    T = zonotope(interval([-.15; -.1],[.15; .1]));

    for t=1:1:horizon
        if mod(t, k1) == 1
            X = zonotope(interval([x(1,1); x(2, 1)],[x(1,1); x(2, 1)]));
            W = get_max_w(X, T, U, A, dt, k1, k2);
        end
        u = stochastic_control(x, max_u, eta);
        supervised_u(3, t) = u;
        if not(isempty(W.vertices))
             [u, super] = supervision(u, W);
        end
        supervised_u(1, t) = u;
        supervised_u(2, t) = super;
        x = move_exact(x, u, dt); 
        supervised_traj(1, t) = x(1, 1);
        supervised_traj(2, t) = x(2, 1);
    end



    
    subplot(3,4, 2*k-1);
    title('position');
    hold on;
    ylim([-1, 1]);
    plot(1:1:horizon, naive_traj(1, :), 'b');
    plot(1:1:horizon, supervised_traj(1, :), 'r');
    plot(1:1:horizon, zeros(1, horizon), 'k--');


    col = 'kr';
    subplot(3, 4, 2*k);
    title('u');
    hold on;
    for t = 1:horizon
       scatter(t, supervised_u(1, t), strcat(col(supervised_u(2, t) + 1),'.'));
       if supervised_u(2, t) == 1
            scatter(t, supervised_u(3, t), 'b.')    
       end
    end
    plot(1:1:horizon, naive_u, 'b');
% plot(1:1:horizon, supervised_u(1, :), 'k');
end






















function [res, super] = supervision(u, w)
    inter = w.vertices;
    sup = inter(2, 2);
    low = inter(2, 1);
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

function res = k_step_forward(x, u, A, dt, k)
    for i = 1:k
        x = A*x + dt*u;
        f = plot(x, [1,2], 'r', 'linewidth', .5);
    end
    res = x;
end
function res = k_step_backward(x, u, A, dt, k)
    for i = 1:k
        x = inv(A) *  (x + -dt*u);
        if mod(i, 2) == 0
            b = plot(x, [1,2], 'b', 'linewidth', 0.5);
        end
    end
    res = x;
end

function w = get_max_w(X, T, U, A, dt, k1, k2)
    % n_points = 1;
    Theta = [0, pi];
    c = zeros(2, length(Theta));
    d = ones(length(Theta), 1);
    for i = 1:length(Theta)
        theta = Theta(i);
        c(1, i) = cos(theta);
        c(2, i) = sin(theta); % n_points direction on the unit circle.
        sumai = [0 0; 0 0];
        l = c(:, i);
        for k = k2:k1+k2-1
            sumai = sumai + A^k;
        end
        invsumai = inv(sumai); % lhs
        % rhs
        rho_s = supportFunc(T, l);
        rho_aiu = 0;
        for k = 0:k2-1
            rho_aiu = rho_aiu + supportFunc(-A^k*U*dt, l);
        end 
        rho_ak1k2s = supportFunc(A^(k1+k2)*X, l);
        rhs = rho_s + rho_aiu - rho_ak1k2s;
        d(i, 1) = rhs;
    end

    w_box = (1/dt)*invsumai*mptPolytope(c', d);
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

function res = pd_control(x, sup_u)
    g = -10.;
    l = 1.;
    m = 1.;
    beta = 1.;
    a = -g/l;
    b = -beta/m;
    u = -(a+1)*x(1,1) - (b+2)*x(2, 1);
    res = clip(u, -sup_u, sup_u);
end

function res = stochastic_control(x, max_u, eta)
    g = -10.;
    l = 1.;
    m = 1.;
    beta = 1.;
    a = -g/l;
    b = -beta/m;
    u = -(a+1)*x(1,1) - (b+2)*x(2, 1);
    res = clip(u + eta * (2 * max_u *rand() - max_u), -max_u, max_u);
end

function res = move_exact(x, u, dt)
    A = [0 1; 10 -1]; 
    inva = inv(A);
    K = x + inva * [0; u];
    x  = expm(dt * A) * K - inva * [0; u]; 
    res = x; 
end

function res = move(x, u, A, dt)
        res = A*x + [0; dt]*u;
end