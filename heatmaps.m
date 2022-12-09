% PARAMS
dt = 0.1;
g = -10.;
l = 1.;
m = 1.;
beta = 1.;
A = [1. dt; -dt*g/l 1-beta/m*dt];
k1 = 17;
k2 = 9;
max_u = 6;
% END PARAMS

close all

for eta = 0
    rng(10*eta);
    h1 = stochastic_heatmap(A, dt, max_u, 20, eta, 25);
    rng(10*eta);
    h2 = supervision_heatmap(A, dt, max_u, k1, k2, 20, eta, 25);
    plot_heatmap(h1, h2);
end

function h = supervision_heatmap(A, dt, max_u, k1, k2, dim, eta, horizon)
    heat = ones(dim);
    X1 = linspace(-1, 1, dim);
    X2 = linspace(-2, 2, dim);
    T = zonotope(interval([0.; 0.],[0.; 0.]));
    U = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
    for i = 1:dim
        for j = 1:dim
            x2 = X2(dim + 1 - i);
            x1 = X1(j);
            x = [x1; x2];
            for t = 1:horizon
                if mod(t, k1) == 1
                    X = zonotope(interval([x(1,1); x(2, 1)],[x(1,1); x(2, 1)]));
                    W = get_max_w(X, T, U, A, dt, k1, k2);
                 end
                 u = stochastic_control(x, max_u, eta);
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
                u = stochastic_control(x, max_u, eta);
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

function plot_heatmap(h1, h2)
    c = get_contour(h1, h2);
    h1c =  h1 + 10 * c;
    figure;
    own = [winter; 1. 1. 1.];
    heatmap(h1c, 'ColorLimits',[0 3], 'Colormap', own); 
end

function res = get_contour(h1, h2)
    mh1 = ones(length(h1));
    mh2 = ones(length(h1));
    for i = 1:length(h1)
        for j = 1:length(h1)
            if h1(i, j) > 1
                mh1(i, j) = 1;
            else
                mh1(i, j) = 0;
            end
            if h2(i, j) > 1
                mh2(i, j) = 1;
            else
                mh2(i, j) = 0;
            end
        end
    end
    res = mh1-mh2;
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

function res = move(x, u, A, dt)
        res = A*x + [0; dt]*u;
end