close all;clear all; format long;
% True values to simulate data
k0 = 1.0; 
k1 = 0.0;
h = 2.0;

% Parameters
N  = 200;
dt = 1e-2;
T  = 3.0;
T_infy = 1.0;
big_U = 400;
L = 1;

% u0fun = @(x) (1/3).*(x) - (1/300);
u0fun = @(x) sin(pi.*x);  % satisfies u(0)=u(1)=0

[x,t, Utrue, full_A, rhs, biot] = forward_heat_varK(k0,k1,N,dt,T,u0fun, T_infy, h, big_U, L);

% M = 10; % Number of sensors
% sensor_idx = round(linspace(2, N, M));
% x_sensors = x(sensor_idx);

% K = 25; % Number of observation times
% obs_times = linspace(0.01, T, K);
% obs_idx = round(obs_times/dt) + 1;
% t_obs = t(obs_idx);

% Y_clean = Utrue(sensor_idx, obs_idx);

% % add noise
% p = 0.02;  % 2% relative noise
% sigma = p * max(abs(Y_clean(:)));

% rng(0); % reproducible
% Y_noisy = Y_clean + sigma * randn(size(Y_clean));

% save
% save('synthetic_data.mat', 'x_sensors','t_obs','Y_noisy','Y_clean', ...
%      'k0_true','k1_true','h_true','sigma','T_infy');

%% --- Plotting a movie --- 
writerObj = VideoWriter('my_animation.mp4', 'MPEG-4');
open(writerObj);
figure;
h = plot(x, Utrue(:,1), 'r.', 'LineWidth', 2);
% pause;
xlabel('x'); ylabel('u(x,t)');
grid on;
ymin = min(Utrue(:));
ymax = max(Utrue(:));
xlim([min(x) max(x)]);
ylim([ymin ymax]);

nt = size(Utrue,2);
% Plot a movie
for i = 1:nt 
    set(h, 'YData', Utrue(:,i));
    xlabel('x'); ylabel('u(x,T)');
    title(sprintf('Variable-k heat equation at t=%.4f', t(i)));
    grid on;
    frame = getframe(gcf); 
    writeVideo(writerObj, frame);
end

close(writerObj);

xi = linspace(0,1,1000);
yi = (T_infy/big_U).*((biot(1).*(xi-1))/(1-biot(1))-1);
plot(xi,yi, 'b--', 'LineWidth',2); hold on;
plot(x, Utrue(:,i))