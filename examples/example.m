addpath('../')

% Create a rank-deficient dataset
N = 3;
n = 12;
m = [6, 5, 4];
A = zeros(sum(m),n);
for i=1:N
    [Ui, Si, Vi]= svd(randn(m(i),n));
    Si(1,1) = 0;
    if i == 1; V1 = Vi; else; Vi = V1; end
    A(1+sum(m(1:i-1)):sum(m(1:i-1))+m(i),:) = Ui*Si*Vi';
end

% Call the HO-GSVD (which calls the HO-CSD after padding A)
[U, S, V, Tau, ~, taumin, taumax] = hogsvd(A, N, m, n);

% Plot the generalized singular values and the eigenvalues of T
figure;
subplot(N+1,1,1);
plot(1:n, (diag(Tau)-taumin)./(taumax-taumin),...
    'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Eigenvalues of T relative to \tau_{min} and \tau_{max}')
xlim([1 n]);
ylim([0 1]);
for i = 1 : N
    subplot(N+1, 1, i+1);
    Si = S(1+(i-1)*n : i*n, :);
    plot(1:n,diag(Si),...
        'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
    xlim([1 n]);
    ylim([0 1]);
    title(sprintf('HO-GSVs of A_%d',i));
end
