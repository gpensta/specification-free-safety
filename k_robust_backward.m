function res = k_robust_backward(X, U, W, A, B, k)
    for i = 1:k
        d = minkDiff(X, W);
        if isempty(d) 
           break 
        else
            X = inv(A) *  (d + -B*U);
        end
    end
    res = X;
end
