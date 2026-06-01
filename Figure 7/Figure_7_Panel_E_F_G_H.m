clear; clc;

F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k_base = [0.000322195 4000 22235 7152.095014 157337.8291 1.702385143 1.835658901 ...
          2.254209123 2.023631011 1.370890126 1.606015384 1.84551856 0.623690942 ...
          0.032968342 0.030786527 0.101382613 0.074767546 0.89477407 0.004803181 ...
          8.167464174 8565.435226 2.814544497 0.055457107 5.621638845 18.85233629 ...
          1.768328496 0.004710273 5.820217585 80527.43275 2.50692042 0.028545884];

mu = 0.0003468419;

tspan = [0, 1e6];
y0 = zeros(9,1);
N = 500;

sRNA_transcription = 1;
gfp_vals = linspace(0, 2*sRNA_transcription, N);

odeOpts = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-10, ...
    'NonNegative', 1:9 ...
);

sp_ann_fixedR = [0, sRNA_transcription/4, sRNA_transcription/2, sRNA_transcription*3/4, 2*sRNA_transcription];
sp_ann_RC     = [0, sRNA_transcription/4, sRNA_transcription/2, sRNA_transcription*3/4, 2*sRNA_transcription];

sp_pred_RC     = [0, sRNA_transcription/2, sRNA_transcription*4/5, sRNA_transcription, 2*sRNA_transcription];
sp_pred_fixedR = [0, sRNA_transcription*0.6, sRNA_transcription*3/4, sRNA_transcription*0.9, 2*sRNA_transcription];

Y_pred_RC = run_sweep(k_base, F, mu, sRNA_transcription, sp_pred_RC, ...
                      gfp_vals, tspan, y0, odeOpts, false, false);

Y_pred_fixedR = run_sweep(k_base, F, mu, sRNA_transcription, sp_pred_fixedR, ...
                          gfp_vals, tspan, y0, odeOpts, true, false);

Y_ann_RC = run_sweep(k_base, F, mu, sRNA_transcription, sp_ann_RC, ...
                     gfp_vals, tspan, y0, odeOpts, false, true);

Y_ann_fixedR = run_sweep(k_base, F, mu, sRNA_transcription, sp_ann_fixedR, ...
                         gfp_vals, tspan, y0, odeOpts, true, true);

Y_pred_RC     = normalise_by_curve_max(Y_pred_RC);
Y_pred_fixedR = normalise_by_curve_max(Y_pred_fixedR);
Y_ann_RC      = normalise_by_curve_max(Y_ann_RC);
Y_ann_fixedR  = normalise_by_curve_max(Y_ann_fixedR);

plot_four_panel(gfp_vals, sRNA_transcription, ...
                Y_pred_RC, sp_pred_RC, ...
                Y_pred_fixedR, sp_pred_fixedR, ...
                Y_ann_RC, sp_ann_RC, ...
                Y_ann_fixedR, sp_ann_fixedR);


function plot_four_panel(gfp_vals, sRNA_t, ...
                         Y_pred_RC, sp_pred_RC, ...
                         Y_pred_fixedR, sp_pred_fixedR, ...
                         Y_ann_RC, sp_ann_RC, ...
                         Y_ann_fixedR, sp_ann_fixedR)

    figure('Color','w');
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    colors = lines(numel(sp_pred_RC));
    for m = 1:numel(sp_pred_RC)
        plot(gfp_vals, Y_pred_RC(:,m), '-', 'LineWidth', 2, 'Color', colors(m,:));
    end
    xline(sRNA_t, '--k', 'LineWidth', 1, 'HandleVisibility','off');
    title('Predation');
    ylabel('Normalised GFPmut3', 'Interpreter','latex');
    format_axes(gca);
    legend(make_legend(sp_pred_RC), 'Interpreter','latex', 'Location','best');
    hold off;

    nexttile; hold on;
    colors = lines(numel(sp_pred_fixedR));
    for m = 1:numel(sp_pred_fixedR)
        plot(gfp_vals, Y_pred_fixedR(:,m), '-', 'LineWidth', 2, 'Color', colors(m,:));
    end
    xline(sRNA_t, '--k', 'LineWidth', 1, 'HandleVisibility','off');
    title('Predation Fixed R');
    format_axes(gca);
    legend(make_legend(sp_pred_fixedR), 'Interpreter','latex', 'Location','best');
    hold off;

    nexttile; hold on;
    colors = lines(numel(sp_ann_RC));
    for m = 1:numel(sp_ann_RC)
        plot(gfp_vals, Y_ann_RC(:,m), '-', 'LineWidth', 2, 'Color', colors(m,:));
    end
    xline(sRNA_t, '--k', 'LineWidth', 1, 'HandleVisibility','off');
    title('Annihilation');
    xlabel('GFPmut3 transcription rate (nM/s)', 'Interpreter','latex');
    ylabel('Normalised GFPmut3', 'Interpreter','latex');
    format_axes(gca);
    legend(make_legend(sp_ann_RC), 'Interpreter','latex', 'Location','best');
    hold off;

    nexttile; hold on;
    colors = lines(numel(sp_ann_fixedR));
    for m = 1:numel(sp_ann_fixedR)
        plot(gfp_vals, Y_ann_fixedR(:,m), '-', 'LineWidth', 2, 'Color', colors(m,:));
    end
    xline(sRNA_t, '--k', 'LineWidth', 1, 'HandleVisibility','off');
    title('Annihilation Fixed R');
    xlabel('GFPmut3 transcription rate (nM/s)', 'Interpreter','latex');
    format_axes(gca);
    legend(make_legend(sp_ann_fixedR), 'Interpreter','latex', 'Location','best');
    hold off;
end


function leg = make_legend(spRNA_t_vals)

    leg = cell(1, numel(spRNA_t_vals));

    for i = 1:numel(spRNA_t_vals)

        if mod(spRNA_t_vals(i), 1) == 0
            spRNA_str = sprintf('%.0f', spRNA_t_vals(i));
        else
            spRNA_str = sprintf('%.3f', spRNA_t_vals(i));
            spRNA_str = regexprep(spRNA_str, '0+$', '');
            spRNA_str = regexprep(spRNA_str, '\.$', '');
        end

        leg{i} = sprintf('$spRNA_t = %s\\ \\mathrm{nM/s}$', spRNA_str);
    end
end


function Ynorm = normalise_by_curve_max(Y)

    Ynorm = Y;

    for j = 1:size(Y,2)

        max_val = max(Y(:,j));

        if max_val <= 0 || ~isfinite(max_val)
            Ynorm(:,j) = 0;
        else
            Ynorm(:,j) = Y(:,j) ./ max_val;
        end
    end
end


function ss_sum_all = run_sweep(k_base, F, mu, sRNA_t, spRNA_t_vals, gfp_vals, tspan, y0, odeOpts, fixedR, annihilation)

    k = k_base;

    if annihilation
        k(18) = 0;
    end

    ss_sum_all = NaN(numel(gfp_vals), numel(spRNA_t_vals));

    for m = 1:numel(spRNA_t_vals)

        spRNA_t = spRNA_t_vals(m);
        yIC = y0;

        for idx = 1:numel(gfp_vals)

            G = gfp_vals(idx);

            rhs = @(t,y) odesys_rc_predation(t, y, k, F, sRNA_t, spRNA_t, G, mu, fixedR);

            [~, Y] = ode15s(rhs, tspan, yIC, odeOpts);

            y_ss = Y(end,:)';
            y_ss = max(y_ss, 0);

            ss_sum_all(idx,m) = max(0, y_ss(2) + y_ss(3));

            yIC = y_ss;
        end
    end
end


function dy = odesys_rc_predation(~, y, k, F, sRNA_t, spRNA_t, G, mu, fixedR)

    if fixedR
        denom = 1;
    else
        denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
        denom = max(denom, 1e-12);
    end

    term1 = (k(1)*y(1)*y(4)*k(2)) / (k(4)*denom);
    term2 = (k(13)*y(5)*y(4)*k(2)) / (k(5)*denom);

    dy = zeros(9,1);

    dy(1) = G - (mu + F(1))*y(1) - term1;
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


function format_axes(ax)

    grid(ax, 'off');
    box(ax, 'off');

    ax.YAxis.Exponent = 0;

    ylim(ax, [0 1.05]);
    yticks(ax, 0:0.25:1);
    ytickformat(ax, '%.2f');
end