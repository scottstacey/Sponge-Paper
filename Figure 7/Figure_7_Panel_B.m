clear; clc;

%% Parameters 
% F = Fixed Parameters
% F = [1 = delta_m, 2 = tau_GFP, 3 = delta_s, 4 = tau_mscarlet] 
F = [0.001699, 0.002818, 0.0003610, 0.0004495];

%% Fitted Parameters
k_base = [1.47E-06	4000	22235	7152.095014	157337.8291	1.702385143	1.835658901	2.254209123	2.023631011	1.370890126	1.606015384	1.84551856	0.623690942	0.032968342	0.030786527	0.101382613	0.074767546	0.89477407	0.004803181	8.167464174	8565.435226	2.814544497	0.055457107	5.621638845	18.85233629	1.768328496	0.004710273	5.820217585	80527.43275	2.50692042	0.028545884];

sRNA_transcription = 4;  
mu                 = 0.0003468419;

% ODE time span
tspan = [0, 1e6];  % seconds

y0 = [0, 0, 0, 0];

% GFP sweep
N        = 500;
gfp_vals = linspace(0, 2*sRNA_transcription, N);

% k1 multipliers/gamma_sm multiplier
k1_mults = [1, 1.5, 2, 5, 10, 100];
colors   = lines(numel(k1_mults));

ss_sum_all = NaN(N, numel(k1_mults));

%% Loop over k(1) variants
for m = 1:numel(k1_mults)
    k = k_base;
    k(1) = k_base(1) * k1_mults(m);
    
    % steady-state for each GFPmut3_transcription
    ss_sum = NaN(size(gfp_vals));
    for idx = 1:numel(gfp_vals)
        Gval = gfp_vals(idx);
        [~, Y] = ode15s(@(t,y) sRNA_ODEs(t,y,k,F, sRNA_transcription, Gval, mu), tspan, y0);
        y_ss   = Y(end,:);
        ss_sum(idx) = y_ss(2) + y_ss(3);
    end
    
    ss_sum_all(:,m) = ss_sum;
end

%% Plot all curves
figure; hold on;

legend_entries = cell(1, numel(k1_mults)+1);

for m = 1:numel(k1_mults)
    plot(gfp_vals, ss_sum_all(:,m), '-', 'LineWidth', 2, 'Color', colors(m,:));
    legend_entries{m} = sprintf('\\gamma_{sm} = %.4g \\gamma_{sm}', k1_mults(m));
end

% dashed vertical line at sRNA_transcription
xline(sRNA_transcription, '--k', 'LineWidth', 1);
legend_entries{end} = 'sRNA transcription';

xlabel('GFPmut3 transcription rate (nM/s)', 'Interpreter','latex');
ylabel('GFPmut3 Steady State (nM)',         'Interpreter','latex');
title('GFPmut3-sRNA Response -- Predation Simple Mass Action', 'Interpreter','latex');

grid off;
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.0f');

legend(legend_entries, 'Location','best');

hold off;

%% ODE function
function dy = sRNA_ODEs(~, y, k, F, sRNA_trans, GFPmut3_trans, mu)
    G      = GFPmut3_trans;
    sRNA_t = sRNA_trans;

    dy = zeros(4,1);

    dy(1) = G - (mu + F(1))*y(1) - k(1)*y(1)*y(4);
    dy(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dy(3) = F(2)*y(2) - mu * y(3);
    dy(4) = sRNA_t - (mu + F(3))*y(4) - 0*k(1)*y(1)*y(4);
end