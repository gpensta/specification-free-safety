clear;

m = 1.;  g = -9.81; l = 1.; beta = 1 ; max_u  = 2; dt = .01; 
U = max_u * interval(-1,1); Uz = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));
X = .001 * zonotope(interval([-1.; -1.],[1. ; 1.]));

% a = (m * g * l * cos(0)) / (m * l * l);
% b = l  / (m * l * l);
% 
% A = [1,    dt;
%     -dt * a,  1 - dt * beta / (m * l * l)]; 
% B = [0;
%      dt * b];
[A, B] = get_lin_model();
n = 40;
Bw = 0;

for i = 1:n
    Bw = Bw + A^(-i) * Uz * dt;
end

S = Bw;
temp = S.interval(); x1 = temp(1); x2 = temp(2);
W = find_w(S, Uz);
T = 0; 

for i = 1:n
    T = T + A^(i-1)*W;
end

Us = minkDiff(Bw, W);
Us =0.85 * (Us & Uz);
F = Us + W;
reB = k_robust_backward(T, Uz, W, A, B, n);

% Plots

close all; figure; hold on
% xlim([-.1, .1]);
% ylim([-.3, .3]);
plot(zonotope(cartProd(x1, x2)), [1,2], 'r', 'linewidth', 2.,'FaceColor', 'none');
plot(W, [1,2], 'b', 'linewidth', 2.,'FaceColor', 'none');
plot(Bw, [1,2], 'g', 'linewidth', 2.,'FaceColor', 'none');
plot(F, [1,2], 'm', 'linewidth', 2.,'FaceColor', 'none');
plot(W, [1,2], 'c', 'linewidth', 2.,'FaceColor', 'none');
% plot(B1, [1,2], 'r', 'linewidth', 2.,'FaceColor', 'none');
plot(T, [1,2], 'b', 'linewidth', 2.,'FaceColor', 'none');
plot(reB, [1,2], 'k', 'linewidth', 0.5,'FaceColor', 'none');







