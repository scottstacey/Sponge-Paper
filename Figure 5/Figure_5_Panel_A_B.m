clear; clc; close all;

F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k = [0.000322195 4000 22235 7152.095014 157337.8291 1.702385143 1.835658901 ...
     2.254209123 2.023631011 1.370890126 1.606015384 1.84551856 0.623690942 ...
     0.032968342 0.030786527 0.101382613 0.074767546 0.89477407 0.004803181 ...
     8.167464174 8565.435226 2.814544497 0.055457107 5.621638845 18.85233629 ...
     1.768328496 0.004710273 5.820217585 80527.43275 2.50692042 0.028545884];

mu = 0.0003468419;

t1 = (11.46 - 10)*3600;
t2 = (18.43 - 10)*3600;
t3 = (24.34 - 10)*3600;
tEnd = (40 - 10)*3600;

Van_1 = 0;
Van_2 = 5.5e3;

OC6_1 = 0;
OC6_2 = 50;

IPTG_1 = 0;
IPTG_vals = [25e3, 100e3, 160e3, 200e3, 1e6];

hexColors = {'#006400', '#228B22', '#32CD32', '#7CFC00', '#ADFF2F'};
colors = cellfun(@hex2rgb, hexColors, 'UniformOutput', false);
colors = vertcat(colors{:});

ode_opts = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-10 ...
);

plot_response( ...
    zeros(9,1), ...
    @(t,y,Van,OC6,IPTG) pSS_02_006_ODEs(t, y, k, F, mu, Van, OC6, IPTG), ...
    [t1, t2, t3, tEnd], ...
    Van_1, Van_2, OC6_1, OC6_2, IPTG_1, IPTG_vals, colors, ode_opts, ...
    'pSS-02-006: Tuning GFPmut3 Through Varying Sponge Induction' ...
);

plot_response( ...
    zeros(7,1), ...
    @(t,y,Van,OC6,IPTG) pSS_02_005_ODEs(t, y, k, F, mu, Van, OC6, IPTG), ...
    [t1, t2, t3, tEnd], ...
    Van_1, Van_2, OC6_1, OC6_2, IPTG_1, IPTG_vals, colors, ode_opts, ...
    'pSS-02-005: Tuning GFPmut3 Through Varying Sponge Induction' ...
);


function plot_response(y0, ode_fun, tSeg, Van_1, Van_2, OC6_1, OC6_2, IPTG_1, IPTG_vals, colors, ode_opts, fig_title)

    t1 = tSeg(1);
    t2 = tSeg(2);
    t3 = tSeg(3);
    tEnd = tSeg(4);

    t_shared = [];
    y_shared = [];

    [t, Y] = ode15s(@(t,y) ode_fun(t, y, Van_1, OC6_1, IPTG_1), [0 t1], y0, ode_opts);
    t_shared = [t_shared; t];
    y_shared = [y_shared; Y];

    [t, Y] = ode15s(@(t,y) ode_fun(t, y, Van_2, OC6_1, IPTG_1), [t1 t2], Y(end,:)', ode_opts);
    t_shared = [t_shared; t(2:end)];
    y_shared = [y_shared; Y(2:end,:)];

    [t, Y] = ode15s(@(t,y) ode_fun(t, y, Van_2, OC6_2, IPTG_1), [t2 t3], Y(end,:)', ode_opts);
    t_shared = [t_shared; t(2:end)];
    y_shared = [y_shared; Y(2:end,:)];

    figure('Color','w'); hold on;

    plot(t_shared/3600, y_shared(:,3), ...
        'LineWidth', 2, ...
        'Color', colors(end,:), ...
        'HandleVisibility', 'off');

    y0_branch = y_shared(end,:)';
    legend_entries = cell(1, numel(IPTG_vals));

    for i = 1:numel(IPTG_vals)

        IPTG = IPTG_vals(i);

        [t_branch, Y_branch] = ode15s( ...
            @(t,y) ode_fun(t, y, Van_2, OC6_2, IPTG), ...
            [t3 tEnd], ...
            y0_branch, ...
            ode_opts ...
        );

        if i == numel(IPTG_vals)
            t_plot = [t3; t_branch(2:end)];
            y_plot = [y0_branch(3); Y_branch(2:end,3)];
        else
            t_plot = t_branch;
            y_plot = Y_branch(:,3);
        end

        plot(t_plot/3600, y_plot, 'LineWidth', 2, 'Color', colors(i,:));

        legend_entries{i} = sprintf('IPTG = %.0f µM', IPTG/1000);
    end

    xlabel('Time (hours)');
    ylabel('Mature GFPmut3 Concentration per Cell (nM)');
    title(fig_title);
    legend(legend_entries, 'Location', 'best');

    grid off; box off;

    ax = gca;
    ax.YAxis.Exponent = 0;
    ytickformat('%.0f');

    hold off;
end


function dY = pSS_02_006_ODEs(~, y, k, F, mu, Van, OC6, IPTG)

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
    denom = max(denom, 1e-12);

    GFP_t = k(19) + (k(20)-k(19))*Van^k(22)/(Van^k(22) + k(21)^k(22));
    sRNA_t = k(23) + (k(24)-k(23))*OC6^k(26)/(OC6^k(26) + k(25)^k(26));
    spRNA_t = k(27) + (k(28)-k(27))*IPTG^k(30)/(IPTG^k(30) + k(29)^k(30));

    mRNA_term = k(1)*y(1)*y(4)*k(2)/(k(4)*denom);
    spRNA_term = k(13)*y(5)*y(4)*k(2)/(k(5)*denom);

    dY = zeros(9,1);

    dY(1) = GFP_t - (mu + F(1))*y(1) - mRNA_term;
    dY(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dY(3) = F(2)*y(2) - mu*y(3);

    dY(4) = sRNA_t ...
          - (mu + F(3))*y(4) ...
          - mRNA_term ...
          - spRNA_term ...
          + k(18)*k(17)*y(6);

    dY(5) = spRNA_t - (mu + k(14))*y(5) - spRNA_term;

    dY(6) = mRNA_term - (k(17) + mu)*y(6);
    dY(7) = spRNA_term - (k(17) + mu)*y(7);

    dY(8) = k(15)*y(5) + k(15)*y(7) - mu*y(8) - F(4)*y(8);
    dY(9) = F(4)*y(8) - mu*y(9);
end


function dY = pSS_02_005_ODEs(~, y, k, F, mu, Van, OC6, IPTG)

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
    denom = max(denom, 1e-12);

    GFP_t = k(19) + (k(20)-k(19))*Van^k(22)/(Van^k(22) + k(21)^k(22));
    sRNA_t = k(23) + (k(24)-k(23))*OC6^k(26)/(OC6^k(26) + k(25)^k(26));
    spRNA_t = k(27) + (k(28)-k(27))*IPTG^k(30)/(IPTG^k(30) + k(29)^k(30));

    mRNA_term = k(1)*y(1)*y(4)*k(2)/(k(4)*denom);
    spRNA_term = k(13)*y(5)*y(4)*k(2)/(k(5)*denom);

    dY = zeros(7,1);

    dY(1) = GFP_t - (mu + F(1))*y(1) - mRNA_term;
    dY(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dY(3) = F(2)*y(2) - mu*y(3);

    dY(4) = sRNA_t ...
          - (mu + F(3))*y(4) ...
          - mRNA_term ...
          - spRNA_term ...
          + k(18)*k(17)*y(6);

    dY(5) = spRNA_t - (mu + k(16))*y(5) - spRNA_term;

    dY(6) = mRNA_term - (k(17) + mu)*y(6);
    dY(7) = spRNA_term - (k(17) + mu)*y(7);
end


function rgb = hex2rgb(hex)

    rgb = sscanf(hex(2:end), '%2x%2x%2x', [3 1])'/255;
end