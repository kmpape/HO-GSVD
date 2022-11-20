clear all
close all
clc

% Create a rank-deficient dataset
N = 3;
n = 12;
m = [6, 5, 4];
A = zeros(sum(m),n);
for i=1:N
    Ai = randn(m(i),n);
    [Ui, Si, Vi]= svd(Ai);
    Si(1,1) = 0;
    if i == 1; V1 = Vi; else; Vi = V1; end
    A(1+sum(m(1:i-1)):sum(m(1:i-1))+m(i),:) = Ui*Si*Vi';
end
% Ai share the same right singular vectors, and the first singular value is
% zero (at least). So dim(range(A))<= 5.

r = rank(A);
rd = n - r;
[UA,SA,VA] = svd(A);
Upad = blkdiag(UA, eye(rd));
Spad = [SA; zeros(rd,r), eye(rd)];
Apad = Upad*Spad*VA';
[Qall, Rall] = qr(Apad, 0);
Nall = N + 1;
mall = [m, rd];
ppi = 1/Nall;

Q=Qall; N=Nall; m=mall; ppi=0.1;
ZEROTOL = 1e-14; EPS_REL_ISO = 1e-3; EPS_SVD_ISO = 1e-3;
warn_eps_ortho_iso = 0.01; warn_eps_eig_iso = 0.01;

[Uall, Sall, Zall, Tauall, Tall] = hocsd(Qall, mall, 'ppi', ppi);
[Emin, Emax] = magic_eigenvalues_T(N, ppi);


% Call the HO-GSVD (which calls the HO-CSD after padding A)
[U, S, V, Tau, taumin, taumax, iso_classes] = hogsvd(Apad, m);
[U_rd, S_rd, V_rd, Tau_rd, taumin_rd, taumax_rd, iso_classes_rd] = hogsvd(A, m(1:end-1));

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
