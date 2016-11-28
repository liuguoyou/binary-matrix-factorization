function [u, s, v] = ssvd(X, lambdaU, lambdaV)

normX = norm(X);
X = X / normX;

[u,~,v] = svds(X,1);

N = 100;

SST = sum(X(:).^2);

for i = 1:N
    uOld = u;
    vOld = v;

    xu = X' * u;
    v = sign(xu) .* max(0, abs(xu) - lambdaV);
    v = v / norm(v);

    xv = X * v;
    u = sign(xv) .* max(0, abs(xv) - lambdaU);
    u = u / norm(u);
    
    if norm(uOld-u, 1) < 1e-4 || norm(vOld-v, 1) < 1e-4
        break
    end
end

s = normX * u' * X * v;

end
