k0 = 1.0; 
k1 = 0.9;
N  = 200;
dt = 1e-4;
T  = 0.10;

u0fun = @(x) sin(pi*x);  % satisfies u(0)=u(1)=0

[x,t,U] = forward_heat_varK(k0,k1,N,dt,T,u0fun);
nt = size(U,2);

writerObj = VideoWriter('my_animation.mp4', 'MPEG-4');
open(writerObj);
figure;
h = plot(x, U(:,1), 'r.', 'LineWidth', 2);
xlabel('x'); ylabel('u(x,t)');
grid on;
ymin = min(U(:));
ymax = max(U(:));
xlim([min(x) max(x)]);
ylim([ymin ymax]);

% Plot a movie
for i = 1:nt
    set(h, 'YData', U(:,i));
    xlabel('x'); ylabel('u(x,T)');
    title(sprintf('Variable-k heat equation at t=%.4f', t(i)));
    grid on;
    frame = getframe(gcf); 
    writeVideo(writerObj, frame);
end

close(writerObj);
