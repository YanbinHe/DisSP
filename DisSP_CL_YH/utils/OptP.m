function P = OptP(n,neighbours)
cvx_begin
cvx_precision default

variable P(n, n) symmetric
W_bar = zeros(n,n);
for i = 1:n
    for j = neighbours{i}
        W_bar = W_bar + P(i,j)/n * GenerateWij(n,i,j);
    end
end

% cost
minimize( lambda_sum_largest(W_bar, 2) );

subject to
for i = 1:n
    pi = 0;
    for j = neighbours{i}
        P(i,j) >= 0;
        pi = pi + P(i,j);
    end
    pi == 1;
end
cvx_end
P = full(P);
end

function Wij = GenerateWij(n,i,j)
ei = zeros(n,1); ei(i) = 1;
ej = zeros(n,1); ej(j) = 1;
Wij = eye(n,n) - 0.5 .* (ei - ej) * (ei - ej)';
end