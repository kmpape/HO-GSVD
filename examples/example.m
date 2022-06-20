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
[~, S, ~, ~, ~, ~, Tau, ~, taumin, taumax, mpad] = hogsvd(A, N, m, n);

if length(mpad)>length(m), is_padded=1; else; is_padded=0; end

% Plot the generalized singular values and the eigenvalues of T
figure;
subplot(N+1+is_padded,1,1);
plot(1:n, (diag(Tau)-taumin)./(taumax-taumin),...
    'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Eigenvalues of T relative to \tau_{min} and \tau_{max}')
xticks(1:1:n); xlim([1 n]); ylim([0 1]);
for i = 1 : N+is_padded
    subplot(N+1+is_padded, 1, i+1);
    Si = S(1+(i-1)*n : i*n, :);
    plot(1:n,diag(Si),...
        'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
    xticks(1:1:n);
    xlim([1 n]);
    ylim([0 1]);
    if is_padded && (i==N+is_padded)
        title(sprintf('HO-GSVs of A_%d (padding)',i));
    else
        title(sprintf('HO-GSVs of A_%d',i));
    end
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
[~, S, V, ~, ~, ~, Tau, ~, taumin, taumax, mpad] = hogsvd(A, N, m, n, 'ACCELERATE', true);

if length(mpad) > length(m)
    N = N+1;
end

% Plot the generalized singular values and the eigenvalues of T
figure;
subplot(N+1,1,1);
plot(1:n, (diag(Tau)-taumin)./(taumax-taumin),...
    'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
title('Eigenvalues of T relative to \tau_{min} and \tau_{max}')
xticks(1:1:n); xlim([1 n]); ylim([0 1]);
for i = 1 : N
    subplot(N+1, 1, i+1);
    Si = S(1+(i-1)*n : i*n, :);
    plot(1:n,diag(Si),...
        'k--','Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k');
    xticks(1:1:n); xlim([1 n]); ylim([0 1]); title(sprintf('HO-GSVs of A_%d',i));
end
