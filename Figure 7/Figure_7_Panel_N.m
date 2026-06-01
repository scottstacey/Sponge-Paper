clear; clc; close all;

data = readtable('Data/BindingRates_RC.xlsx');

k1_variants  = data.("k_1_");
k13_variants = data.("k_13_");

nVar = numel(k1_variants);

mu = 0.0003468419;

t_half  = 7*3600;
t_final = 14*3600;

GFP_tx = 0.5;

F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k_base = [0.000322195 4000 22235 7152.095014 157337.8291 1.702385143 1.835658901 ...
          2.254209123 2.023631011 1.370890126 1.606015384 1.84551856 0.623690942 ...
          0.032968342 0.030786527 0.101382613 0.074767546 0.89477407 0.004803181 ...
          8.167464174 8565.435226 2.814544497 0.055457107 5.621638845 18.85233629 ...
          1.768328496 0.004710273 5.820217585 80527.43275 2.50692042 0.028545884];

y0 = zeros(9,1);

nGrid = 300;
maxTX = 2*GFP_tx;

sRNA_vals  = linspace(0, maxTX, nGrid);
spRNA_vals = linspace(0, maxTX, nGrid);

odeOpts = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-10, ...
    'NonNegative', 1:9 ...
);

pool = gcp('nocreate');

if isempty(pool)
    parpool;
end

nPts = nGrid*nGrid;
meanChangeVec = NaN(nPts, 1);

parfor linIdx = 1:nPts

    [is, ip] = ind2sub([nGrid, nGrid], linIdx);

    sRNA_tx  = sRNA_vals(is);
    spRNA_tx = spRNA_vals(ip);

    rhs1 = @(t,y) odeFun_directTX(t, y, k_base, F, GFP_tx, sRNA_tx, spRNA_tx, mu);

    [~, Y1] = ode15s(rhs1, [0 t_half], y0, odeOpts);

    y_switch = max(Y1(end,:)', 0);
    Y_switch = y_switch(2) + y_switch(3);

    if Y_switch <= 0
        continue;
    end

    changes = NaN(nVar, 1);

    for v = 1:nVar

        k_mod = k_base;
        k_mod(1) = k1_variants(v);
        k_mod(13) = k13_variants(v);

        rhs2 = @(t,y) odeFun_directTX(t, y, k_mod, F, GFP_tx, sRNA_tx, spRNA_tx, mu);

        [~, Y2] = ode15s(rhs2, [t_half t_final], y_switch, odeOpts);

        y_end = max(Y2(end,:)', 0);
        Y_end = y_end(2) + y_end(3);

        changes(v) = (Y_end/Y_switch) - 1;
    end

    meanChangeVec(linIdx) = mean(changes, 'omitnan');
end

meanChange = reshape(meanChangeVec, nGrid, nGrid);

figure('Color','w');

imagesc(spRNA_vals, sRNA_vals, meanChange);
set(gca, 'YDir', 'normal');

colorbar;

xlabel('spRNA transcription rate (nM/s)');
ylabel('sRNA transcription rate (nM/s)');
title(sprintf('Mean normalised expression change after mutation, GFP TX = %.2f nM/s', GFP_tx));

grid off; box off;


function dy = odeFun_directTX(~, y, k, F, GFP_t, sRNA_t, spRNA_t, mu)

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
    denom = max(denom, 1e-12);

    mRNA_term = k(1)*y(1)*y(4)*k(2)/(k(4)*denom);
    spRNA_term = k(13)*y(5)*y(4)*k(2)/(k(5)*denom);

    dy = zeros(9,1);

    dy(1) = GFP_t - (mu + F(1))*y(1) - mRNA_term;
    dy(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dy(3) = F(2)*y(2) - mu*y(3);

    dy(4) = sRNA_t ...
          - (mu + F(3))*y(4) ...
          - mRNA_term ...
          - spRNA_term ...
          + k(18)*k(17)*y(6);

    dy(5) = spRNA_t - (mu + k(14))*y(5) - spRNA_term;

    dy(6) = mRNA_term - (k(17) + mu)*y(6);
    dy(7) = spRNA_term - (k(17) + mu)*y(7);

    dy(8) = k(15)*y(5) + k(15)*y(7) - mu*y(8) - F(4)*y(8);
    dy(9) = F(4)*y(8) - mu*y(9);
end