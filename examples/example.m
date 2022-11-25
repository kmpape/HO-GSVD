addpath('../')

%% Example 1: Dataset with Rank-Deficient A
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

% Call the HO-GSVD, which calls the HO-CSD after padding A
[U, S, V, Tau, taumin, taumax, iso_classes] = hogsvd(A, m);

% Plot the generalized singular values and the eigenvalues of T
figure;
subplot(N+1,1,1);
plot(1:n, (diag(Tau)-taumin)./(taumax-taumin),...
    'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Eigenvalues \tau_k of T_\pi as (\tau_k-\tau_{min})/(\tau_{max}-\tau_{min})')
xticks(1:1:n); xlim([1 n]); ylim([0 1]);
for i = 1 : N
    subplot(N+1, 1,i+1);
    Si = S(1+(i-1)*n : i*n, :);    
    plot(1:n,diag(Si),...
        'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
    title(sprintf('Generalized Singular Values \\sigma_{%d,k} of A_{%d}',i,i))
    xticks(1:1:n);
    xlim([1 n]);
    ylim([0 1]);
end

%% Example 2: Full-Rank Dataset A with common and isolated subspace
N = 3;
n = 4;
m = [4, 3, 3];
A = zeros(sum(m),n);
for i=1:N
    [Ui, Si, Vi]= svd(randn(m(i),n));
    Si(1,1) = 1;
    if i == 1; V1 = Vi; else; Vi = V1; end
    A(1+sum(m(1:i-1)):sum(m(1:i-1))+m(i),:) = Ui*Si*Vi';
end

% Call the HO-GSVD (which calls the HO-CSD after padding A)
[U, S, V, Tau, taumin, taumax, iso_classes] = hogsvd(A, m);

% Plot the generalized singular values and the eigenvalues of T
figure;
subplot(N+1,1,1);
plot(1:n, (diag(Tau)-taumin)./(taumax-taumin),...
    'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Eigenvalues \tau_k of T_\pi as (\tau_k-\tau_{min})/(\tau_{max}-\tau_{min})')
xticks(1:1:n); xlim([1 n]); ylim([0 1]);
for i = 1 : N
    subplot(N+1, 1, i+1);
    Si = S(1+(i-1)*n : i*n, :);
    plot(1:n,diag(Si),...
        'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
    title(sprintf('Generalized Singular Values \\sigma_{%d,k} of A_{%d}',i,i))
    xticks(1:1:n); xlim([1 n]); ylim([0 1]);
end
