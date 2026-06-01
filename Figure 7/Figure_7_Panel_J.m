clear; clc;

%% Constants
t1 = 2 * 3600; 
t2 = 4 * 3600; 
y0 = zeros(9,1);

Van_init  = [0.5, 0.5, 0.5, 0.5];
Van_later = [0.5, 0.5, 1.1915, 0.5];

OC6_init  = [0.1, 0.1, 0.1, 0.1];
OC6_later = [0.1, 0,   0.1, 0];

IPTG_init  = [0, 0, 0, 0];
IPTG_later = [4, 0, 0, 4];

labels = {'Sponge Induced', 'sRNA Switched Off', 'GFPmut3 Induced More', 'sRNA Switched Off + Sponge Induced'};
colors = lines(4);

% Parameters
F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k_base = [0.000322195	4000	22235	7152.095014	157337.8291	1.702385143	1.835658901 ...
          2.254209123	2.023631011	1.370890126	1.606015384	1.84551856	0.623690942 ...
          0.032968342	0.030786527	0.101382613	0.074767546	0.89477407	0.004803181 ...
          8.167464174	8565.435226	2.814544497	0.055457107	5.621638845	18.85233629 ...
          1.768328496	0.004710273	5.820217585	80527.43275	2.50692042	0.028545884];

mu = 0.0003468419;

fixedR  = false;
inducer = false;

%% Time vectors
tvec1 = linspace(0,  t1,     300);
tvec2 = linspace(t1, t1+t2,  300);
t_all_sec = [tvec1, tvec2(2:end)];
t_all_hr  = t_all_sec / 3600;

%% Rise time metrics setup
rt = table('Size',[4 3], ...
    'VariableTypes', {'string','double','double'}, ...
    'VariableNames', {'Condition','t50_hr','t90_hr'});
rt.Condition = string(labels(:));

%% Solve and Plot
figure('Color','w'); hold on;

for i = 1:4
    % First segment
    [~, Y1] = ode45(@(t,y) odeFun(t,y,k_base,F,Van_init(i),OC6_init(i),IPTG_init(i),mu,fixedR,inducer), tvec1, y0);

    % Second segment
    [~, Y2] = ode45(@(t,y) odeFun(t,y,k_base,F,Van_later(i),OC6_later(i),IPTG_later(i),mu,fixedR,inducer), tvec2, Y1(end,:)');

    % Combine
    Y_full = [Y1; Y2(2:end,:)];
    y_out  = Y_full(:,3);

    % Plot
    plot(t_all_hr, y_out, 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', labels{i});

    % 50% and 90% rise times after induction
    t50_abs = transition_time_frac(t_all_sec, y_out, t1, 0.50);
    t90_abs = transition_time_frac(t_all_sec, y_out, t1, 0.90);

    rt.t50_hr(i) = (t50_abs - t1) / 3600;
    rt.t90_hr(i) = (t90_abs - t1) / 3600;
end


xline(t1/3600, '--k', 'LineWidth', 1.5, 'DisplayName', 'Inducer Switch');
xlabel('Time (hours)');
ylabel('GFPmut3 per Cell (nM)');
title('Response Time: Step Change in GFPmut3');
legend('Location','best');
grid off; box off;
ax = gca; ax.YAxis.Exponent = 0;
ytickformat('%.0f');
hold off;

disp(rt);


function t_hit = transition_time_frac(t_sec, y, t_switch_sec, frac)

    idx = find(t_sec >= t_switch_sec, 1, 'first');
    if isempty(idx) || idx == numel(t_sec)
        t_hit = NaN; return;
    end

    y0 = y(idx);
    yF = y(end);
    dy = yF - y0;

    if ~isfinite(y0) || ~isfinite(yF) || abs(dy) < 1e-12
        t_hit = NaN; return;
    end

    target = y0 + frac * dy;

    if dy > 0
        j = find(y(idx:end) >= target, 1, 'first');
    else
        j = find(y(idx:end) <= target, 1, 'first');
    end
    if isempty(j)
        t_hit = NaN; return;
    end
    j = j + idx - 1;

    if j == idx
        t_hit = t_sec(j); return;
    end

    t1 = t_sec(j-1); t2 = t_sec(j);
    y1 = y(j-1);     y2 = y(j);

    if ~isfinite(y1) || ~isfinite(y2) || abs(y2 - y1) < 1e-15
        t_hit = t_sec(j); return;
    end

    alpha = (target - y1) / (y2 - y1);
    alpha = min(max(alpha, 0), 1);
    t_hit = t1 + alpha * (t2 - t1);
end

%% ODE system
function dy = odeFun(~, y, k, F, Van, OC6, IPTG, mu, fixedR, inducer)

    if inducer
        sRNA_t  = (k(23) + (k(24)-k(23))*((OC6^k(26))/(OC6^k(26) + k(25)^k(26))));
        GFP_t   = (k(19) + (k(20)-k(19))*((Van^k(22))/(Van^k(22) + k(21)^k(22))));
        spRNA_t = (k(27) + (k(28)-k(27))*((IPTG^k(30))/(IPTG^k(30) + k(29)^k(30))));
    else
        sRNA_t  = OC6;
        GFP_t   = Van;
        spRNA_t = IPTG;
    end

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