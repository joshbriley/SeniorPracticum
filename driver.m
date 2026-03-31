close all;clear all; format long;
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');

% True values to simulate data
alpha = 1.0;
alpha_a = 0.0; 
alpha_b = 2.0;
n_MC = 5000;

% Uniform discrete prior over sampled alpha candidates
prior = ones(n_MC, 1) ./ n_MC;
alpha_dist = alpha_a + (alpha_b - alpha_a).*rand(n_MC, 1);
log_likelihood = zeros(n_MC, 1);
h = 20.0;

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

M = 30; % Number of sensors
sensor_idx = round(linspace(2, N, M));
x_sensors = x(sensor_idx);

K = 30; % Number of observation times
obs_times = linspace(0.01, T, K);
obs_idx = round(obs_times/dt) + 1;
t_obs = t(obs_idx);

Y_clean = Utrue(sensor_idx, obs_idx);

% add noise
p = 0.1; 
sigma = p * max(abs(Y_clean(:)));

% rng(0); % reproducible
Y_noisy = Y_clean + sigma * randn(size(Y_clean));

% Store predicted sensor values for each candidate alpha so the posterior
% can be updated one observation time at a time.
Y_model_all = zeros(n_MC, M, K);
for i = 1:n_MC
    alpha_i = alpha_dist(i);

    % Solve forward problem with candidate alpha
    [~,~, U_i] = forward_heat_varK(alpha_i, N, dt, T, ...
        u0fun, T_infy, h, big_U, L);

    % Extract model predictions at sensor locations
    Y_model_all(i,:,:) = U_i(sensor_idx, obs_idx);
end

% Sequentially update the posterior as each observation time arrives.
alpha_map_history = zeros(K, 1);
alpha_mean_history = zeros(K, 1);
post_by_time = zeros(n_MC, K);

for k = 1:K
    misfit_k = bsxfun(@minus, Y_noisy(:,k)', Y_model_all(:,:,k));
    log_likelihood = log_likelihood ...
        - (1/(2*sigma^2)) * sum(misfit_k.^2, 2);

    log_post_k = log(prior) + log_likelihood;
    log_post_k = log_post_k - max(log_post_k);
    post_k = exp(log_post_k);
    post_by_time(:,k) = post_k;

    weights_k = post_k ./ sum(post_k);
    [~, idx_map_k] = max(weights_k);
    alpha_map_history(k) = alpha_dist(idx_map_k);
    alpha_mean_history(k) = sum(weights_k .* alpha_dist);
end

post = post_by_time(:, end);

% Sort by alpha for numerical integration
[alpha_sorted, sort_idx] = sort(alpha_dist);
post_sorted = post(sort_idx);

% Numerical integral via trapezoid rule
Z = trapz(alpha_sorted, post_sorted);  % normalizing constant
post_normalized = post_sorted / Z;     
post_normalized = post_normalized;

% MAP estimate
[~, idx_max] = max(post_normalized);
alpha_MAP = alpha_sorted(idx_max);
post_max = post_normalized(idx_max);
alpha_mean = sum((post ./ sum(post)) .* alpha_dist);

disp(['True alpha: ', num2str(alpha)])
disp(['Alpha MAP estimate: ', num2str(alpha_MAP)])
disp(['Alpha posterior mean: ', num2str(alpha_mean)])
disp(['Max PDF value: ', num2str(post_max)])

save('synthetic_data.mat', 'x_sensors','t_obs','Y_noisy','Y_clean', ...
    'alpha','h','sigma','T_infy','alpha_map_history','alpha_mean_history');

figure;
scatter(alpha_sorted, post_normalized, 30, 'b','filled');
ylim([0 1]);
xlabel('$\alpha$', 'Fontsize', 20);
ylabel('Posterior Probability', 'Fontsize', 20);
fig_title = "Posterior distribution of $\alpha$ with relative noise = " +  p;
title(fig_title);
grid on; hold on;
xline(alpha, '-r', 'LineWidth', 2)
legend('Probability Distribution Function', 'True value of $\alpha$', 'Fontsize', 12)
set(gca, 'FontSize', 14);
filename = "figs/posterior_scatter_" + p + ".png"; 
exportgraphics(gcf, filename, 'Resolution', 750);

figure;
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
exportgraphics(gcf, 'figs/alpha_estimate_over_time.png', 'Resolution', 750);

%% --- Plotting a movie --- (Comment for speed) 
% writerObj = VideoWriter('my_animation.mp4', 'MPEG-4');
% open(writerObj);
% figure;
% h = plot(x, Utrue(:,1), 'r.', 'LineWidth', 2);
% xlabel('x'); ylabel('u(x,t)');
% grid on;
% ymin = min(Utrue(:));
% ymax = max(Utrue(:));
% xlim([min(x) max(x)]);
% ylim([ymin ymax]);

% nt = size(Utrue,2);
% % Plot a movie
% for i = 1:nt 
%     set(h, 'YData', Utrue(:,i));
%     xlabel('x'); ylabel('u(x,T)');
%     title(sprintf('Heat Equation at t=%.4f', t(i)));
%     grid on;
%     frame = getframe(gcf); 
%     writeVideo(writerObj, frame);
%     pause(1);
%     filename = sprintf('figs/forward_heat_%.1f.png', i);
%     exportgraphics(gcf, filename, 'Resolution', 500);
% end

% close(writerObj);
% hold on;
% xi = linspace(0,1,1000);
% yi = (T_infy/big_U).*((biot.*(xi-1))/(1-biot)-1);
% plot(xi,yi, 'b--', 'LineWidth',2); hold on;
% legend('Approx', 'True Steady State')
