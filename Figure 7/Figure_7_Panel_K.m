clear; clc; close all;

F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k_base = [0.000322195 4000 22235 7152.095014 157337.8291 1.702385143 1.835658901 ...
          2.254209123 2.023631011 1.370890126 1.606015384 1.84551856 0.623690942 ...
          0.032968342 0.030786527 0.101382613 0.074767546 0.89477407 0.004803181 ...
          8.167464174 8565.435226 2.814544497 0.055457107 5.621638845 18.85233629 ...
          1.768328496 0.004710273 5.820217585 80527.43275 2.50692042 0.028545884];

mu = 0.0003468419;
fixedR = false;

tPert_hr = 6;
tEnd_hr  = 18;

tPert = tPert_hr * 3600;
tEnd  = tEnd_hr  * 3600;

tgrid = linspace(0, tEnd, 2001);
tgrid = sort(unique([tgrid, tPert]));
time_hr = tgrid / 3600;

idxPert = find(tgrid == tPert, 1, 'first');

y0 = zeros(9,1);

odeOpts = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-10, ...
    'NonNegative', 1:9 ...
);

GFP0  = 0.5;
sRNA0 = 0.05;
sp0   = 0.05;

rng(1);

noiseSD      = 0.1;
smoothWin_hr = 1;

u0 = sp0;

u = u0 * ones(size(tgrid));
postMask = tgrid >= tPert;
tPost = tgrid(postMask);

dt_hr = median(diff(tgrid)) / 3600;
winN = max(3, round(smoothWin_hr / dt_hr));

z = randn(size(tPost));
z = smoothdata(z, 'movmean', winN);
z = z - mean(z);

uPost = u0 + noiseSD * z;
uPost = max(0, uPost);

u(postMask) = uPost;

GFP_A   = GFP0  * ones(size(tgrid));
sRNA_A  = sRNA0 * ones(size(tgrid));
spRNA_A = sp0   * ones(size(tgrid));
spRNA_A(postMask) = u(postMask);

GFP_B   = GFP0  * ones(size(tgrid));
sRNA_B  = sRNA0 * ones(size(tgrid));
spRNA_B = zeros(size(tgrid));
sRNA_B(postMask) = u(postMask);

YA = simulate_rc(tgrid, y0, k_base, F, mu, fixedR, odeOpts, GFP_A, sRNA_A, spRNA_A);
YB = simulate_rc(tgrid, y0, k_base, F, mu, fixedR, odeOpts, GFP_B, sRNA_B, spRNA_B);

GFP_A = YA(:,3);
GFP_B = YB(:,3);

GFP_A = GFP_A ./ max(GFP_A(idxPert), 1e-12);
GFP_B = GFP_B ./ max(GFP_B(idxPert), 1e-12);

figure('Color','w');
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');

nexttile([3 1]); hold on;
plot(time_hr, GFP_A, 'LineWidth', 2, 'DisplayName','Noisy spRNA input');
plot(time_hr, GFP_B, 'LineWidth', 2, 'DisplayName','Noisy sRNA input');
xline(tPert_hr, '--k', 'LineWidth', 1.5, 'HandleVisibility','off');
ylabel('Normalised GFPmut3');
title('Noisy input filtering');
legend('Location','best');
grid off; box off;
ax = gca;
ax.YAxis.Exponent = 0;

nexttile; hold on;
plot(time_hr, u, 'LineWidth', 2, 'DisplayName','Noisy input');
yline(u0, '--', 'LineWidth', 1.2, 'DisplayName','Mean input');
xline(tPert_hr, '--k', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('Time (hours)');
ylabel('Input TX (nM/s)');
legend('Location','best');
grid off; box off;
ax = gca;
ax.YAxis.Exponent = 0;


function Y = simulate_rc(tgrid, y0, k, F, mu, fixedR, odeOpts, GFP_t, sRNA_t, spRNA_t)

    FG = griddedInterpolant(tgrid(:), GFP_t(:),   'linear', 'nearest');
    FS = griddedInterpolant(tgrid(:), sRNA_t(:),  'linear', 'nearest');
    FP = griddedInterpolant(tgrid(:), spRNA_t(:), 'linear', 'nearest');

    rhs = @(t,y) odeFun_core(t, y, k, F, mu, fixedR, FG, FS, FP);

    [~, Y] = ode15s(rhs, tgrid, y0, odeOpts);
end


function dy = odeFun_core(t, y, k, F, mu, fixedR, FG, FS, FP)

    GFP_t   = FG(t);
    sRNA_t  = FS(t);
    spRNA_t = FP(t);

    if fixedR
        denom = 1;
    else
        denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
        denom = max(denom, 1e-12);
    end

    term1 = (k(1)*y(1)*y(4)*k(2)) / (k(4)*denom);
    term2 = (k(13)*y(5)*y(4)*k(2)) / (k(5)*denom);

    dy = zeros(9,1);

    dy(1) = GFP_t - (mu + F(1))*y(1) - term1;
    dy(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dy(3) = F(2)*y(2) - mu*y(3);

    dy(4) = sRNA_t ...
          - (mu + F(3))*y(4) ...
          - term1 ...
          - term2 ...
          + k(18)*k(17)*y(6);

    dy(5) = spRNA_t ...
          - (mu + k(14))*y(5) ...
          - term2;

    dy(6) = term1 - (k(17)+mu)*y(6);
    dy(7) = term2 - (k(17)+mu)*y(7);

    dy(8) = k(15)*y(5) + k(15)*y(7) - mu*y(8) - F(4)*y(8);
    dy(9) = F(4)*y(8) - mu*y(9);
end