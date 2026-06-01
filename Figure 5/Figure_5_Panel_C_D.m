function Figure_5_Panel_C_D()

clear; clc; close all;

F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k = [0.000322195 4000 22235 7152.095014 157337.8291 1.702385143 1.835658901 ...
     2.254209123 2.023631011 1.370890126 1.606015384 1.84551856 0.623690942 ...
     0.032968342 0.030786527 0.101382613 0.074767546 0.89477407 0.004803181 ...
     8.167464174 8565.435226 2.814544497 0.055457107 5.621638845 18.85233629 ...
     1.768328496 0.004710273 5.820217585 80527.43275 2.50692042 0.028545884];

mu = 0.0003468419;

OC6_vec_nM = linspace(0, 120, 100);
IPTG_vec_nM = linspace(0, 300000, 200);
IPTG_vec_uM = IPTG_vec_nM/1000;

OC6_vec_nM_line = linspace(0, 120, 400);
txIPTG_line_vec = [0, 5];

ss.T0 = 1e5;
ss.maxExt = 6;
ss.tol_dydt = 1e-8;
ss.maxStep = 50;

odeOpts = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-10, ...
    'MaxStep', ss.maxStep ...
);

output_fun = @(y) y(2) + y(3);
eps0 = 1e-12;

[Y0, yss0] = steady_state_output( ...
    @(t,y) rhs_circ2_inducer(t, y, 0, 0, mu, F, k), ...
    zeros(9,1), ...
    output_fun, ...
    ss, ...
    odeOpts ...
);

Y0 = max(Y0, eps0);

Y2 = nan(numel(IPTG_vec_nM), numel(OC6_vec_nM));

yRowStart = yss0;

for iy = 1:numel(IPTG_vec_nM)

    IPTG = IPTG_vec_nM(iy);
    yIC = yRowStart;

    for ix = 1:numel(OC6_vec_nM)

        OC6 = OC6_vec_nM(ix);

        [Yij, yssij] = steady_state_output( ...
            @(t,y) rhs_circ2_inducer(t, y, OC6, IPTG, mu, F, k), ...
            yIC, ...
            output_fun, ...
            ss, ...
            odeOpts ...
        );

        Y2(iy, ix) = Yij;
        yIC = yssij;
    end

    yRowStart = yIC;
end

FC2 = Y2/Y0;

climInd = prctile(FC2(isfinite(FC2)), [2 98]);
climInd(1) = max(0, climInd(1));

figure(2); clf;
set(gcf, 'Color', 'w');

imagesc(OC6_vec_nM, IPTG_vec_uM, FC2);
set(gca, 'YDir', 'normal');

xlabel('OC6 (nM)');
ylabel('IPTG (\muM)');
title('Tuneability of spRNA regulated gene expression');

caxis(climInd);

cb = colorbar;
cb.Label.String = 'Relative GFPmut3 Expression (1 = no sRNA/spRNA induction)';

hold on;
contour(OC6_vec_nM, IPTG_vec_uM, FC2, 10, 'k');

FC2_mixed_lines = nan(numel(txIPTG_line_vec), numel(OC6_vec_nM_line));

for iline = 1:numel(txIPTG_line_vec)

    txIPTG = txIPTG_line_vec(iline);

    [Y0_line, yIC] = steady_state_output( ...
        @(t,y) rhs_circ2_mixed(t, y, 0, txIPTG, mu, F, k), ...
        zeros(9,1), ...
        output_fun, ...
        ss, ...
        odeOpts ...
    );

    Y0_line = max(Y0_line, eps0);

    for ix = 1:numel(OC6_vec_nM_line)

        OC6 = OC6_vec_nM_line(ix);

        [Yij, yssij] = steady_state_output( ...
            @(t,y) rhs_circ2_mixed(t, y, OC6, txIPTG, mu, F, k), ...
            yIC, ...
            output_fun, ...
            ss, ...
            odeOpts ...
        );

        FC2_mixed_lines(iline, ix) = Yij/Y0_line;
        yIC = yssij;
    end
end

figure(5); clf;
set(gcf, 'Color', 'w');
hold on;

for iline = 1:numel(txIPTG_line_vec)

    plot( ...
        OC6_vec_nM_line, ...
        FC2_mixed_lines(iline,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('spRNA tx = %.1f', txIPTG_line_vec(iline)) ...
    );
end

xlabel('OC6 (nM)');
ylabel('Relative steady-state output (normalised to OC6 = 0)');
title('Dose response with fixed direct spRNA transcription rates');
legend('Location', 'best');
grid on;

end


function dydt = rhs_circ2_inducer(~, y, OC6, IPTG, mu, F, k)

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
    denom = max(denom, 1e-12);

    term1 = k(1)*y(1)*y(4)*k(2)/(k(4)*denom);
    term2 = k(13)*y(5)*y(4)*k(2)/(k(5)*denom);

    prodOC6 = k(23) + (k(24)-k(23))*OC6^k(26)/(OC6^k(26) + k(25)^k(26));
    prodIPTG = k(27) + (k(28)-k(27))*IPTG^k(30)/(IPTG^k(30) + k(29)^k(30));

    dydt = zeros(9,1);

    dydt(1) = 1 - (mu + F(1))*y(1) - term1;
    dydt(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dydt(3) = F(2)*y(2) - mu*y(3);

    dydt(4) = prodOC6 ...
            - (mu + F(3))*y(4) ...
            - term1 ...
            - term2 ...
            + k(18)*k(17)*y(6);

    dydt(5) = prodIPTG - (mu + k(14))*y(5) - term2;

    dydt(6) = term1 - (k(17) + mu)*y(6);
    dydt(7) = term2 - (k(17) + mu)*y(7);

    dydt(8) = k(15)*y(5) + k(15)*y(7) - mu*y(8) - F(4)*y(8);
    dydt(9) = F(4)*y(8) - mu*y(9);
end


function dydt = rhs_circ2_mixed(~, y, OC6, txIPTG, mu, F, k)

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
    denom = max(denom, 1e-12);

    term1 = k(1)*y(1)*y(4)*k(2)/(k(4)*denom);
    term2 = k(13)*y(5)*y(4)*k(2)/(k(5)*denom);

    prodOC6 = k(23) + (k(24)-k(23))*OC6^k(26)/(OC6^k(26) + k(25)^k(26));
    prodIPTG = txIPTG;

    dydt = zeros(9,1);

    dydt(1) = 1 - (mu + F(1))*y(1) - term1;
    dydt(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dydt(3) = F(2)*y(2) - mu*y(3);

    dydt(4) = prodOC6 ...
            - (mu + F(3))*y(4) ...
            - term1 ...
            - term2 ...
            + k(18)*k(17)*y(6);

    dydt(5) = prodIPTG - (mu + k(14))*y(5) - term2;

    dydt(6) = term1 - (k(17) + mu)*y(6);
    dydt(7) = term2 - (k(17) + mu)*y(7);

    dydt(8) = k(15)*y(5) + k(15)*y(7) - mu*y(8) - F(4)*y(8);
    dydt(9) = F(4)*y(8) - mu*y(9);
end


function [Yss, yss] = steady_state_output(rhs, y0, outfun, ss, odeOpts)

    T = ss.T0;
    yIC = y0;

    for ext = 0:ss.maxExt

        [~, y] = ode15s(rhs, [0 T], yIC, odeOpts);
        yss = y(end,:)';

        if norm(rhs(T, yss), 2) <= ss.tol_dydt
            Yss = outfun(yss);
            return;
        end

        yIC = yss;
        T = 2*T;
    end

    Yss = outfun(yss);
end