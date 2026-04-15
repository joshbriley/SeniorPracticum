close all;clear all; format long;
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');

if ~exist('figs', 'dir')
    mkdir('figs');
end

% True values to simulate data
alpha = 1.0;
alpha_a = 1e-3; 
alpha_b = 2.0;
h = 20.0;
h_a = 5.0;
h_b = 35.0;
n_MC = 7500;

% Uniform discrete prior over sampled (alpha, h) candidates
prior = ones(n_MC, 1) ./ n_MC;
alpha_dist = alpha_a + (alpha_b - alpha_a).*rand(n_MC, 1);
h_dist = h_a + (h_b - h_a).*rand(n_MC, 1);
log_likelihood = zeros(n_MC, 1);

% Parameters
N  = 300;
dt = 1e-2;
T  = 1.0;
T_infy = 1.0;
big_U = 400;
L = 1;

% u0fun = @(x) (1/3).*(x) - (1/300);
u0fun = @(x) cos(pi.*x);  

[x,t, Utrue, full_A, rhs, biot] = forward_heat_varK(alpha, N, dt, T, u0fun, T_infy, h, big_U, L);

M = 10; % Number of sensors
sensor_idx = round(linspace(2, N, M));
x_sensors = x(sensor_idx);

K = 25; % Number of observation times
obs_times = linspace(0.01, T, K);
obs_idx = round(obs_times/dt) + 1;
t_obs = t(obs_idx);

Y_clean = Utrue(sensor_idx, obs_idx);

% add noise
p = 0.1; 
sigma = p * max(abs(Y_clean(:)));

% rng(0); % reproducible
Y_noisy = Y_clean + sigma * randn(size(Y_clean));

% Store predicted sensor values for each candidate parameter pair so the posterior
% can be updated one observation time at a time.
Y_model_all = nan(n_MC, M, K);
for i = 1:n_MC
    alpha_i = alpha_dist(i);
    h_i = h_dist(i);

    % Solve forward problem with candidate (alpha, h)
    [~,~, U_i] = forward_heat_varK(alpha_i, N, dt, T, ...
        u0fun, T_infy, h_i, big_U, L);

    if all(isfinite(U_i(:)))
        % Extract model predictions at sensor locations
        Y_model_all(i,:,:) = U_i(sensor_idx, obs_idx);
    end
end

% Sequentially update the posterior as each observation time arrives.
alpha_map_history = zeros(K, 1);
alpha_mean_history = zeros(K, 1);
h_map_history = zeros(K, 1);
h_mean_history = zeros(K, 1);
post_by_time = zeros(n_MC, K);

for k = 1:K
    pred_k = squeeze(Y_model_all(:,:,k));
    misfit_k = bsxfun(@minus, pred_k, Y_noisy(:,k)');
    valid_k = all(isfinite(misfit_k), 2);

    log_likelihood(~valid_k) = -Inf;
    if any(valid_k)
        log_likelihood(valid_k) = log_likelihood(valid_k) ...
            - (1/(2*sigma^2)) * sum(misfit_k(valid_k,:).^2, 2);
    end

    log_post_k = log(prior) + log_likelihood;
    finite_mask_k = isfinite(log_post_k);
    if ~any(finite_mask_k)
        error('All parameter candidates became invalid by observation step %d.', k);
    end

    post_k = zeros(n_MC, 1);
    max_log_post_k = max(log_post_k(finite_mask_k));
    post_k(finite_mask_k) = exp(log_post_k(finite_mask_k) - max_log_post_k);
    post_by_time(:,k) = post_k;

    sum_post_k = sum(post_k);
    if ~isfinite(sum_post_k) || sum_post_k <= 0
        error('Posterior normalization failed at observation step %d.', k);
    end
    weights_k = post_k ./ sum_post_k;
    [~, idx_map_k] = max(weights_k);
    alpha_map_history(k) = alpha_dist(idx_map_k);
    alpha_mean_history(k) = sum(weights_k .* alpha_dist);
    h_map_history(k) = h_dist(idx_map_k);
    h_mean_history(k) = sum(weights_k .* h_dist);
end

post = post_by_time(:, end);
sum_post = sum(post);
if ~isfinite(sum_post) || sum_post <= 0
    error('Final posterior normalization failed.');
end
weights = post ./ sum_post;

% MAP estimate
[post_max, idx_max] = max(weights);
alpha_MAP = alpha_dist(idx_max);
h_MAP = h_dist(idx_max);
alpha_mean = sum(weights .* alpha_dist);
h_mean = sum(weights .* h_dist);

disp(['True alpha: ', num2str(alpha)])
disp(['True h: ', num2str(h)])
disp(['Alpha MAP estimate: ', num2str(alpha_MAP)])
disp(['h MAP estimate: ', num2str(h_MAP)])
disp(['Alpha posterior mean: ', num2str(alpha_mean)])
disp(['h posterior mean: ', num2str(h_mean)])
disp(['Max posterior weight: ', num2str(post_max)])

save('synthetic_data.mat', 'x_sensors','t_obs','Y_noisy','Y_clean', ...
    'alpha','h','sigma','T_infy','alpha_map_history','alpha_mean_history', ...
    'h_map_history','h_mean_history','alpha_dist','h_dist','weights');

figure;
scatter(alpha_dist, h_dist, 30, weights, 'filled');
xlabel('$\alpha$', 'Fontsize', 20);
ylabel('$h$', 'Fontsize', 20);
fig_title = sprintf('Posterior samples for $(\\alpha,h)$ with relative noise = %.3g', p);
title(fig_title);
grid on; hold on;
plot(alpha, h, 'rp', 'MarkerSize', 14, 'LineWidth', 2);
plot(alpha_MAP, h_MAP, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
cb = colorbar;
cb.Label.String = 'Posterior weight';
legend('Posterior samples', 'True value', 'MAP', 'Fontsize', 12, ...
    'Location', 'best');
set(gca, 'FontSize', 14);
filename = sprintf('figs/posterior_scatter_%.3g.png', p);
exportgraphics(gcf, filename, 'Resolution', 750);

figure;
subplot(2,1,1);
alpha_edges = linspace(alpha_a, alpha_b, 41);
alpha_centers = 0.5 * (alpha_edges(1:end-1) + alpha_edges(2:end));
alpha_bin_width = alpha_edges(2) - alpha_edges(1);
alpha_prior_counts = histc(alpha_dist, alpha_edges);
alpha_prior_counts(end-1) = alpha_prior_counts(end-1) + alpha_prior_counts(end);
alpha_prior_counts = alpha_prior_counts(1:end-1);
alpha_prior_density = alpha_prior_counts ./ (sum(alpha_prior_counts) * alpha_bin_width);
alpha_post_counts = weighted_histc(alpha_dist, weights, alpha_edges);
alpha_post_density = alpha_post_counts ./ alpha_bin_width;
stairs(alpha_centers, alpha_prior_density, 'LineWidth', 1.5); hold on;
bar(alpha_centers, alpha_post_density, 1.0, 'FaceAlpha', 0.35, ...
    'EdgeColor', 'none');
xline(alpha, ':', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2]);
xlabel('$\alpha$', 'Fontsize', 20);
ylabel('Density', 'Fontsize', 20);
title('Prior and posterior marginals for $\alpha$');
legend('Prior samples', 'Posterior marginal', 'True value', ...
    'Fontsize', 12, 'Location', 'best');
grid on;

subplot(2,1,2);
h_edges = linspace(h_a, h_b, 41);
h_centers = 0.5 * (h_edges(1:end-1) + h_edges(2:end));
h_bin_width = h_edges(2) - h_edges(1);
h_prior_counts = histc(h_dist, h_edges);
h_prior_counts(end-1) = h_prior_counts(end-1) + h_prior_counts(end);
h_prior_counts = h_prior_counts(1:end-1);
h_prior_density = h_prior_counts ./ (sum(h_prior_counts) * h_bin_width);
h_post_counts = weighted_histc(h_dist, weights, h_edges);
h_post_density = h_post_counts ./ h_bin_width;
stairs(h_centers, h_prior_density, 'LineWidth', 1.5); hold on;
bar(h_centers, h_post_density, 1.0, 'FaceAlpha', 0.35, ...
    'EdgeColor', 'none');
xline(h, ':', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2]);
xlabel('$h$', 'Fontsize', 20);
ylabel('Density', 'Fontsize', 20);
title('Prior and posterior marginals for $h$');
legend('Prior samples', 'Posterior marginal', 'True value', ...
    'Fontsize', 12, 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);
exportgraphics(gcf, 'figs/prior_posterior_marginals.png', ...
    'Resolution', 750);

figure;
subplot(2,1,1);
plot(t_obs, alpha_map_history, 'b-o', 'LineWidth', 1.5); hold on;
plot(t_obs, alpha_mean_history, 'k--', 'LineWidth', 1.5);
yline(alpha, '-r', 'LineWidth', 2);
xlabel('Observation time', 'Fontsize', 20);
ylabel('$\alpha$ estimate', 'Fontsize', 20);
title('Sequential estimates of $\alpha$');
grid on;
legend('MAP', 'Posterior mean', 'True value', 'Fontsize', 12, ...
    'Location', 'best');
set(gca, 'FontSize', 14);

subplot(2,1,2);
plot(t_obs, h_map_history, 'b-o', 'LineWidth', 1.5); hold on;
plot(t_obs, h_mean_history, 'k--', 'LineWidth', 1.5);
yline(h, '-r', 'LineWidth', 2);
xlabel('Observation time', 'Fontsize', 20);
ylabel('$h$ estimate', 'Fontsize', 20);
title('Sequential estimates of $h$');
grid on;
legend('MAP', 'Posterior mean', 'True value', 'Fontsize', 12, ...
    'Location', 'best');
set(gca, 'FontSize', 14);
exportgraphics(gcf, 'figs/parameter_estimate_over_time.png', 'Resolution', 750);

%% --- Plotting a movie --- (Comment for speed) 
% writerObj = VideoWriter('my_animation.mp4', 'MPEG-4');
% open(writerObj);
% figure;
% set(gca, 'FontSize', 14);
% h_plot = plot(x, Utrue(:,1), 'r.', 'LineWidth', 2);
% xlabel('x', 'Fontsize', 20);
% ylabel('u(x,T)', 'Fontsize', 20);
% grid on;
% ymin = min(Utrue(:));
% ymax = max(Utrue(:));
% xlim([min(x) max(x)]);
% ylim([ymin ymax]);

% nt = size(Utrue,2);
% % Plot a movie
% for i = 1:nt 
%     set(h_plot, 'YData', Utrue(:,i));
%     set(gca, 'FontSize', 14);
%     xlabel('x', 'Fontsize', 16);
%     ylabel('u(x,T)', 'Fontsize', 16);
%     title(sprintf('Heat Equation at t=%.4f', t(i)));
%     grid on;
%     frame = getframe(gcf); 
%     writeVideo(writerObj, frame);
%     pause(1);
%     filename = sprintf('figs/forward_heat_%04d.png', i);
%     exportgraphics(gcf, filename, 'Resolution', 500);
% end

% close(writerObj);
% hold on;
% xi = linspace(0,1,1000);
% yi = (T_infy/big_U).*((biot.*(xi-1))/(1-biot)-1);
% plot(xi,yi, 'b--', 'LineWidth',2); hold on;
% legend('Approx', 'True Steady State')

figure; plot(x, Utrue(:,2), 'r.', 'LineWidth', 2)
xlim([0 1]); ylim([-0.91 0.9])
grid on; hold on;
plot(x_sensors, Y_noisy(:,1), 'bo', 'LineWidth', 2)
xlabel('x'); ylabel('u(x,T)');title('Heat Equation at t=0.0100')
exportgraphics(gcf, 'ForSolvObs_half.png', 'Resolution', 700);
set(gca, 'FontSize', 16);
legend('Forward Solver', 'Synthetic Observational Data', 'Fontsize', 14, ...
    'Location', 'best');

function counts = weighted_histc(samples, sample_weights, edges)
counts = zeros(length(edges) - 1, 1);
for j = 1:length(counts)
    if j < length(counts)
        in_bin = samples >= edges(j) & samples < edges(j+1);
    else
        in_bin = samples >= edges(j) & samples <= edges(j+1);
    end
    counts(j) = sum(sample_weights(in_bin));
end
end
