clear;
max_u  = 2; n = 14;
U = max_u * interval(-1,1); Uz = max_u * zonotope(interval([0.; -1.],[0. ; 1.]));

[A, B] = get_lin_model();
Bw = 0;

for i = 1:n
    Bw = Bw + A^(-i) * B * Uz;
end

S = Bw;
temp = S.interval(); x1 = temp(1); x2 = temp(2);
W =  find_w(S, Uz) + 0 * zonotope(interval([-1.; -1.],[1. ; 1.]));
T = 0; 

for i = 1:n
    T = T + A^(i-1)*W;
end

if in(temp,T)

    diff = minkDiff(Bw, W);

    Us = diff & Uz;

%     Us = .85 * (Us & Uz);
    F = Us + W;
    reB = k_robust_backward(T, Uz, W, A, B, n);
    
    % Plots
    close all; figure; hold on

    plot(Bw, [1,2], 'g', 'linewidth', 2.,'FaceColor', 'none');
    plot(zonotope(cartProd(x1, x2)), [1,2], 'r', 'linewidth', 2.,'FaceColor', 'none');
    plot(W, [1,2], 'b', 'linewidth', 2.,'FaceColor', 'none');
    plot(T, [1,2], 'b', 'linewidth', 2.,'FaceColor', 'none');
%     plot(F, [1,2], 'm', 'linewidth', 2.,'FaceColor', 'none');
    % plot(B1, [1,2], 'r', 'linewidth', 2.,'FaceColor', 'none');
%     plot(reB, [1,2], 'k', 'linewidth', 0.5,'FaceColor', 'none');
    % legend({'$\mathcal{B}_{40}(\mathcal{T}_{num2str(n)})$', '$\mathcal{S}$', '$\mathcal{W}$', '$\mathcal{T}_{40}$', '$\mathcal{F}_1(0)$'},'Interpreter','latex') % R2018a and earlier
    legend({['$\mathcal{B}_{' num2str(n) '}(\mathcal{T}_{' num2str(n) '})$'], '$\mathcal{S}$', '$\mathcal{W}$', ['$\mathcal{T}_{'  num2str(n)  '}$'], '$\mathcal{F}_1(0)$'},'Interpreter','latex') 

end




