clear; clc; close all;

F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k = [0.000322195 4000 22235 7152.095014 157337.8291 1.702385143 1.835658901 ...
     2.254209123 2.023631011 1.370890126 1.606015384 1.84551856 0.623690942 ...
     0.032968342 0.030786527 0.101382613 0.074767546 0.89477407 0.004803181 ...
     8.167464174 8565.435226 2.814544497 0.055457107 5.621638845 18.85233629 ...
     1.768328496 0.004710273 5.820217585 80527.43275 2.50692042 0.028545884];

mu = 0.0003468419;
scale = 10;

tEnd = 12*3600;
tMid = tEnd/2;

tgrid = linspace(0, tEnd, 1201);
tgrid = sort(unique([tgrid, tMid]));
time_hr = tgrid/3600;

y0 = zeros(11,1);

odeOpts = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-10, ...
    'NonNegative', 1:11 ...
);

GFP_t = 0.1;

sRNA_pre = 0;
sRNA_post = 0.1;

spRNA_on = 0.1;
spRNA_off = 0;

rhs_on  = @(t,y) Sponge_1_4_ODEs(t, y, k, F, mu, scale, GFP_t, sRNA_pre, sRNA_post, spRNA_on, tMid);
rhs_off = @(t,y) Sponge_1_4_ODEs(t, y, k, F, mu, scale, GFP_t, sRNA_pre, sRNA_post, spRNA_off, tMid);

[~, Y_on]  = ode15s(rhs_on, tgrid, y0, odeOpts);
[~, Y_off] = ode15s(rhs_off, tgrid, y0, odeOpts);

GFP_on  = Y_on(:,2)  + Y_on(:,3);
GFP_off = Y_off(:,2) + Y_off(:,3);

Off_on  = Y_on(:,9)  + Y_on(:,10);
Off_off = Y_off(:,9) + Y_off(:,10);

figure('Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
plot(time_hr, GFP_on, 'LineWidth', 2, 'DisplayName','+ spRNA');
plot(time_hr, GFP_off, 'LineWidth', 2, 'DisplayName','- spRNA');
xline(tMid/3600, '--k', 'LineWidth', 1.5, 'HandleVisibility','off');
ylabel('GFPmut3 Expression');
title('On-Target Gene Expression Change +/- spRNA');
legend('Location','best');
grid off; box off;

nexttile; hold on;
plot(time_hr, Off_on, 'LineWidth', 2, 'DisplayName','+ spRNA');
plot(time_hr, Off_off, 'LineWidth', 2, 'DisplayName','- spRNA');
xline(tMid/3600, '--k', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('Time (hours)');
ylabel('Off-Target Gene Expression');
title('Off-Target Gene Expression Change +/- spRNA');
legend('Location','best');
grid off; box off;


function dY = Sponge_1_4_ODEs(t, y, k, F, mu, scale, GFP_t, sRNA_pre, sRNA_post, spRNA_t, tMid)

    spRNA_t = spRNA_t*10;

    if t < tMid
        sRNA_t = sRNA_pre;
    else
        sRNA_t = sRNA_post*10;
    end

    OffT_t = GFP_t;
    Off_t_K = k(4)*10;

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5) + (y(8)*y(4))/Off_t_K;
    denom = max(denom, 1e-12);

    term_mRNA = k(1)*k(2)*y(1)*y(4)/(k(4)*denom);
    term_spRNA = k(13)*k(2)*y(4)*y(5)/(k(5)*denom);
    term_offTarget = (k(1)/scale)*k(2)*y(4)*y(8)/(Off_t_K*denom);

    dY = zeros(11,1);

    dY(1)  = GFP_t - (mu + F(1))*y(1) - term_mRNA;
    dY(2)  = k(31)*y(1) - (mu + F(2))*y(2);
    dY(3)  = F(2)*y(2) - mu*y(3);

    dY(4)  = sRNA_t ...
           - (mu + F(3))*y(4) ...
           - term_mRNA ...
           - term_spRNA ...
           - term_offTarget ...
           + k(18)*k(17)*y(6) ...
           + k(18)*k(17)*y(11);

    dY(5)  = spRNA_t - (mu + k(14))*y(5) - term_spRNA;

    dY(6)  = term_mRNA - (k(17) + mu)*y(6);
    dY(7)  = term_spRNA - (k(17) + mu)*y(7);

    dY(8)  = OffT_t - (mu + F(1))*y(8) - term_offTarget;
    dY(9)  = k(31)*y(8) - (mu + F(2))*y(9);
    dY(10) = F(2)*y(9) - mu*y(10);

    dY(11) = term_offTarget - (k(17) + mu)*y(11);
end