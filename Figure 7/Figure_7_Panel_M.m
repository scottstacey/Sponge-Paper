clear; clc; close all;

data = readtable('Data/BindingRates_RC.xlsx');

variant_names = data.("sRNAVariant");
k1_variants   = data.("k_1_");
k13_variants  = data.("k_13_");

F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k_base = [0.000322195 4000 22235 7152.095014 157337.8291 1.702385143 1.835658901 ...
          2.254209123 2.023631011 1.370890126 1.606015384 1.84551856 0.623690942 ...
          0.032968342 0.030786527 0.101382613 0.074767546 0.89477407 0.004803181 ...
          8.167464174 8565.435226 2.814544497 0.055457107 5.621638845 18.85233629 ...
          1.768328496 0.004710273 5.820217585 80527.43275 2.50692042 0.028545884];

mu = 0.0003468419;

t_half  = 7*3600;
t_final = 14*3600;

t1 = linspace(0, t_half, 300);
t2_full = linspace(t_half, t_final, 301);
t2 = t2_full(2:end);

time = [t1, t2];
time_hr = time/3600;

idx_half = find(time == t_half, 1, 'first');

Van1 = 5753; OC61 = 30; IPTG1 = 200000;
Van2 = 6653; OC62 = 30; IPTG2 = 0;

y0 = zeros(9,1);

nVar = numel(k1_variants);
nT = numel(time);

out_with = zeros(nVar, nT);
out_without = zeros(nVar, nT);

for i = 1:nVar

    k_mod = k_base;
    k_mod(1) = k1_variants(i);
    k_mod(13) = k13_variants(i);

    [~, Y1] = ode45(@(t,y) odeFun(t, y, k_base, F, Van1, OC61, IPTG1, mu), t1, y0);
    [~, Y2] = ode45(@(t,y) odeFun(t, y, k_mod, F, Van1, OC61, IPTG1, mu), t2_full, Y1(end,:)');

    Y = [Y1; Y2(2:end,:)];
    out = (Y(:,2) + Y(:,3))';

    out_with(i,:) = out ./ max(out(idx_half), eps);

    [~, Y1] = ode45(@(t,y) odeFun(t, y, k_base, F, Van2, OC62, IPTG2, mu), t1, y0);
    [~, Y2] = ode45(@(t,y) odeFun(t, y, k_mod, F, Van2, OC62, IPTG2, mu), t2_full, Y1(end,:)');

    Y = [Y1; Y2(2:end,:)];
    out = (Y(:,2) + Y(:,3))';

    out_without(i,:) = out ./ max(out(idx_half), eps);
end

mean_with = mean(out_with, 1);
std_with = std(out_with, 0, 1);

mean_without = mean(out_without, 1);
std_without = std(out_without, 0, 1);

figure('Color','w'); hold on;

fill([time_hr fliplr(time_hr)], ...
     [mean_with + std_with fliplr(mean_with - std_with)], ...
     [0.2 0.6 1], ...
     'FaceAlpha', 0.2, ...
     'EdgeColor', 'none', ...
     'HandleVisibility', 'off');

fill([time_hr fliplr(time_hr)], ...
     [mean_without + std_without fliplr(mean_without - std_without)], ...
     [1 0.4 0.4], ...
     'FaceAlpha', 0.2, ...
     'EdgeColor', 'none', ...
     'HandleVisibility', 'off');

plot(time_hr, mean_with, 'b', 'LineWidth', 2, 'DisplayName', 'With Sponge');
plot(time_hr, mean_without, 'r', 'LineWidth', 2, 'DisplayName', 'Without Sponge');

xline(t_half/3600, '--k', 'LineWidth', 1.5, 'DisplayName', 'Mutation');

xlabel('Time (hours)');
ylabel('Mean Normalised GFPmut3 Expression Across Mutants');
title('Effect of Single Nucleotide Mutation to sRNA on Gene Expression +/- spRNA');
legend('Location', 'best');
grid off; box off;


function dy = odeFun(~, y, k, F, Van, OC6, IPTG, mu)

    dy = zeros(9,1);

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
    denom = max(denom, 1e-12);

    GFP_t = k(19) + (k(20)-k(19)) * Van^k(22) / (Van^k(22) + k(21)^k(22));
    sRNA_t = k(23) + (k(24)-k(23)) * OC6^k(26) / (OC6^k(26) + k(25)^k(26));
    spRNA_t = k(27) + (k(28)-k(27)) * IPTG^k(30) / (IPTG^k(30) + k(29)^k(30));

    mRNA_term = k(1)*y(1)*y(4)*k(2) / (k(4)*denom);
    spRNA_term = k(13)*y(5)*y(4)*k(2) / (k(5)*denom);

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