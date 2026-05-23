%% Major Assignment 02: MCMC for Chemical Transformation Kinetics
% Final Implementation: No Toolboxes Required
clc; clear; close all; rng(42);

%% --- 1. KINETIC MODEL ---
% Analytical solution for Parent -> Metabolite -> Sink
function C = model_func(theta, t)
    kp = theta(1); km = theta(2); M0 = theta(3); c = theta(4);
    if abs(kp - km) < 1e-7
        parent = M0 * exp(-kp * t);
        met    = c * M0 * kp * t .* exp(-kp * t);
    else
        parent = M0 * exp(-kp * t);
        met    = c * M0 * kp * (exp(-km*t) - exp(-kp*t)) / (kp - km);
    end
    C = [parent, met];
end

%% --- 2. DATA SETUP ---
t_data   = [0,5,10,20,30,50,70,100,150,200,250,300,365]';
true_A   = [0.058, 0.023, 194, 0.224];   % Scenario A: Non-zero km
true_B   = [0.058, 0.000, 194, 0.224];   % Scenario B: Zero km
sigma_rel = 0.05;

noise_add = @(y) y .* (1 + sigma_rel * randn(size(y)));
C_obs_A   = max(0, noise_add(model_func(true_A, t_data)));
C_obs_B   = max(0, noise_add(model_func(true_B, t_data)));

% [FIGURE 1]: Simulated Data
figure(1);
subplot(1,2,1);
plot(t_data, C_obs_A(:,1), 'ro', 'MarkerSize', 7, 'LineWidth', 1.5); hold on;
plot(t_data, C_obs_A(:,2), 'bs', 'MarkerSize', 7, 'LineWidth', 1.5);
title('Scenario A (km > 0)'); xlabel('Time (d)'); ylabel('Concentration (g/L)');
legend({'Parent','Met'},'Location','northeast'); grid on;

subplot(1,2,2);
plot(t_data, C_obs_B(:,1), 'ro', 'MarkerSize', 7, 'LineWidth', 1.5); hold on;
plot(t_data, C_obs_B(:,2), 'bs', 'MarkerSize', 7, 'LineWidth', 1.5);
title('Scenario B (km = 0)'); xlabel('Time (d)'); ylabel('Concentration (g/L)');
legend({'Parent','Met'},'Location','northeast'); grid on;
saveas(gcf, 'fig1_simulated_data.png');

%% --- 3. LOG-POSTERIOR ---
function lp = log_posterior(vec, t, y)
    kp=vec(1); km=vec(2); M0=vec(3); c=vec(4);
    sp=exp(vec(5)); sm=exp(vec(6));           % sample on log(sigma) scale
    if kp<=0 || km<0 || M0<=0 || c<0 || c>1
        lp = -Inf; return;
    end
    pred = model_func([kp, km, M0, c], t);
    n    = length(t);
    lp   = -0.5*sum(((y(:,1)-pred(:,1))/sp).^2) - n*log(sp) ...
           -0.5*sum(((y(:,2)-pred(:,2))/sm).^2) - n*log(sm);
end

%% --- 4. METROPOLIS-HASTINGS SAMPLER (shared settings) ---
n_iter  = 200000;   % synced with report
burn    = 50000;    % synced with report
thin    = 10;
% Reduce proposal SD to improve acceptance rate
prop_sd = [0.0005, 0.0002, 0.5, 0.002, 0.015, 0.015]; %prop_sd = [0.0015, 0.0008, 1.5, 0.008, 0.04, 0.04];

function [post, acc_rate] = run_mcmc(init, n_iter, burn, thin, prop_sd, t, y)
    samples   = zeros(n_iter, 6);
    samples(1,:) = init;
    cur_lp    = log_posterior(init, t, y);
    acc       = 0;
    for i = 2:n_iter
        prop    = samples(i-1,:) + prop_sd .* randn(1,6);
        prop_lp = log_posterior(prop, t, y);
        if log(rand()) < prop_lp - cur_lp
            samples(i,:) = prop;
            cur_lp       = prop_lp;
            acc          = acc + 1;
        else
            samples(i,:) = samples(i-1,:);
        end
    end
    post     = samples(burn+1:thin:end, :);
    acc_rate = 100 * acc / n_iter;
end

%% --- 5. SCENARIO A: Non-zero km ---
fprintf('Running Scenario A ...\n');
init_A           = [0.05, 0.02, 180, 0.2, log(5), log(5)];
[post_A, accA]   = run_mcmc(init_A, n_iter, burn, thin, prop_sd, t_data, C_obs_A);
fprintf('Scenario A acceptance rate: %.1f%%\n', accA);

% Sort for quantile
sorted_kmA = sort(post_A(:,2));
q5_A = sorted_kmA(max(1, round(0.05*length(sorted_kmA))));
fprintf('Scenario A 5%%-quantile of km: %.4f d^-1\n', q5_A);

% [FIGURE 2]: Posterior histograms — fix ylim issue by drawing line after hist
figure(2);
param_names = {'k_p', 'k_m', 'M_0', 'c'};
true_vals   = true_A;
for i = 1:4
    subplot(2,2,i);
    h = histogram(post_A(:,i), 50, 'FaceColor', [0.2 0.6 0.8], ...
                  'EdgeColor', 'none');
    hold on;
    % Draw true-value line AFTER histogram so ylim is already set
    yl = ylim;
    line([true_vals(i) true_vals(i)], yl, ...
         'Color','r','LineStyle','--','LineWidth',2);
    xlabel(param_names{i}); ylabel('Count');
    title(['Posterior of ', param_names{i}]);
    grid on;
end
sgtitle('Scenario A: Posterior Distributions');
saveas(gcf, 'fig2_posterior_A.png');

% [FIGURE 3]: MCMC Predictions with visible 90% CI
n_pred = 500;
idx    = randperm(size(post_A,1), n_pred);
t_fine = linspace(0, 365, 200)';
all_pred = zeros(200, 2, n_pred);
for k = 1:n_pred
    all_pred(:,:,k) = model_func(post_A(idx(k),1:4), t_fine);
end

% Sort along sample dimension for quantile bands
sorted_pred = sort(all_pred, 3);
lo  = sorted_pred(:, :, round(0.05*n_pred));
hi  = sorted_pred(:, :, round(0.95*n_pred));
mn  = mean(all_pred, 3);

figure(3);
species_name = {'Parent Fit', 'Metabolite Fit'};
for s = 1:2
    subplot(1,2,s);
    % Fill CI band first so the mean line sits on top
    fill([t_fine; flipud(t_fine)], [lo(:,s); flipud(hi(:,s))], ...
         [0.75 0.85 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    plot(t_fine, mn(:,s),    'b-',  'LineWidth', 2);
    plot(t_data, C_obs_A(:,s), 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
    xlabel('Time (d)'); ylabel('Concentration (g/L)');
    title(species_name{s});
    legend({'90% CI', 'Mean', 'Data'}, 'Location', 'northeast');
    grid on;
end
sgtitle('Scenario A: MCMC Predictions');
saveas(gcf, 'fig3_predictions_A.png');

%% --- 6. SCENARIO B: Zero km ---
fprintf('\nRunning Scenario B ...\n');
init_B           = [0.05, 0.001, 180, 0.2, log(5), log(5)];
[post_B, accB]   = run_mcmc(init_B, n_iter, burn, thin, prop_sd, t_data, C_obs_B);
fprintf('Scenario B acceptance rate: %.1f%%\n', accB);

sorted_kmB = sort(post_B(:,2));
q5_B = sorted_kmB(max(1, round(0.05*length(sorted_kmB))));
fprintf('Scenario B 5%%-quantile of km: %.2e d^-1\n', q5_B);

% [FIGURE 4]: Posterior of km (Zero Case)
figure(4);
histogram(post_B(:,2), 60, 'FaceColor', [0.8 0.3 0.3], 'EdgeColor','none');
hold on;
yl = ylim;
line([q5_B q5_B], yl, 'Color','g','LineWidth',2.5);
text(q5_B + 5e-5, yl(2)*0.6, sprintf(' 5%% Q = %.5f', q5_B), 'FontSize', 11);
xlabel('k_m (d^{-1})'); ylabel('Count');
title('Scenario B: Posterior of k_m  (true k_m = 0)');
grid on;
saveas(gcf, 'fig4_posterior_B_km.png');

%% --- 7. FIGURE 5: MCMC vs Normal Approximation ---
% MAP via fmincon (no toolbox-free alternative needed — fmincon is base MATLAB)
neg_lp = @(x) -log_posterior(x, t_data, C_obs_B);
lb = [1e-6, 0,   1,  0, -5, -5];
ub = [1,    1, 500,  1,  5,  5];
opts = optimoptions('fmincon','Display','off','MaxIterations',2000);
x_map = fmincon(neg_lp, init_B, [],[],[],[], lb, ub, [], opts);

% Numerical Hessian (finite differences)
d = 1e-5;
H = zeros(6);
f0 = neg_lp(x_map);
for ii = 1:6
    for jj = 1:6
        ei = zeros(1,6); ei(ii) = 1;
        ej = zeros(1,6); ej(jj) = 1;
        H(ii,jj) = (neg_lp(x_map+d*ei+d*ej) - neg_lp(x_map+d*ei) ...
                   - neg_lp(x_map+d*ej) + f0) / d^2;
    end
end
Sigma   = inv(H);
std_km  = sqrt(abs(Sigma(2,2)));    % abs() guards against tiny negative rounding errors

% Plot range: centre on MAP km, go ±4 sigma but clamp left at -0.01
km_map  = x_map(2);
x_left  = max(-0.015, km_map - 4*std_km);
x_right = km_map + 5*std_km;
x_v     = linspace(x_left, x_right, 300);
norm_pdf = (1/(std_km*sqrt(2*pi))) * exp(-0.5*((x_v - km_map)/std_km).^2);

figure(5);
% Histogram normalised to PDF so both curves share the same y-axis
histogram(post_B(:,2), 80, 'Normalization','pdf', ...
          'FaceColor',[0.5 0.7 0.9], 'EdgeColor','none', 'FaceAlpha', 0.75);
hold on;
plot(x_v, norm_pdf, 'r-', 'LineWidth', 2.5);
% Mark the negative region to emphasise the problem
yl = ylim;
patch([-0.015 0 0 -0.015], [0 0 yl(2) yl(2)], [1 0.8 0.8], ...
      'FaceAlpha',0.3,'EdgeColor','none');
text(-0.012, yl(2)*0.8, 'Physically impossible', ...
     'FontSize', 9, 'Color', [0.7 0 0]);
xline(0, 'k--', 'LineWidth', 1.5);
xlim([-0.015 x_right]);
xlabel('k_m (d^{-1})'); ylabel('Probability Density');
title('Comparison: MCMC Posterior vs Normal Approximation');
legend({'MCMC (true posterior)', 'Normal approximation (misleading)'}, ...
       'Location','northeast');
grid on;
saveas(gcf, 'fig5_comparison_mcmc_vs_normal.png');

%% --- 8. SUMMARY ---
fprintf('\n========== SUMMARY ==========\n');
fprintf('Scenario A acceptance rate : %.1f%%\n', accA);
fprintf('Scenario B acceptance rate : %.1f%%\n', accB);
fprintf('Scenario A 5%%-quantile km : %.4f d^-1  (should be > 0)\n', q5_A);
fprintf('Scenario B 5%%-quantile km : %.2e d^-1  (should be ~ 0)\n', q5_B);
fprintf('Normal approx 5%%-quantile : %.4f d^-1  (negative = impossible)\n', ...
        km_map - 1.645*std_km);
fprintf('All figures saved.\n');