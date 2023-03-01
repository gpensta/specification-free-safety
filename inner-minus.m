% Constants and variables
Z = zonotope([1 2 2 1; 0 0 2 2]);
Z2 = .5*zonotope(interval([-1.; -1.],[1. ; 1.]));
V = 0.5 * [-1 1 1 -1; -1 -1 1 1];
Zd = innerMinkDiff(Z, V);
close all;
figure; hold on;
plot(Z,[1 2], 'b');
plot(Zd,[1 2], 'r');
plot(Zd+Z2,[1 2], 'g');
plot(Z2,[1 2], 'k');

function res = innerMinkDiff(Z, V)
    g = Z.generators;
    N = size(g, 2);
    b = ones(N, 1);
    m = size(V, 2); % Number of vertices in V
    n = size(V, 1); % Dimension of disturbance set W
    c = Z.center; % Initialize c
    dx = N*(m+1)+n; % x = N * m + m + n = m * (N + 1) + n
    % Objective function
    f = zeros([zeros(1, N*m) ones(1, N) zeros(1, 2)]); % b = 1 is the same
    % Constraints
    Aeq = zeros(m*n, dx);
    beq = zeros(m*n, 1);
    A = zeros(N*m, dx);
    B = zeros(2*m*N, 1);
    for j = 1:m
        for i = 1:N
            Aeq(2*j-1, N*(j-1)+i) = g(1, i);
            Aeq(2*j-1, dx-1) = 1;
            Aeq(2*j, N*(j-1)+i) = g(2, i);
            Aeq(2*j, dx) = 1;
        end
        beq(2*j-1) = V(1, j);
        beq(2*j) = V(2, j);
    end

    for i = 1:N
        for j = 1:m
            A((i-1)*m+j, i+(j-1)*N) = 1;
            A((i-1)*m+j, N*m + i) = -1;
        end
    end
    A = [A; -A(1:N*m, 1:N*m) A(:,N*m+1:dx)]; % because it is absolute value on theta
    %%% End Constraints
    [x, fval, exitflag] = linprog(f, A, B, Aeq, beq);
    gr =(1 - x(N*m+1:N*m+N))' .*g;
    cr = c - x(dx-1:dx);
    res = zonotope(cr, gr);