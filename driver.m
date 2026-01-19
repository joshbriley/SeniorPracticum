k0 = 1.0; 
k1 = 0.1;     % must satisfy k0 > |k1|
N  = 200;
dt = 1e-4;
T  = 0.05;

u0fun = @(x) sin(pi*x);  % satisfies u(0)=u(1)=0

[x,t,U] = forward_heat_varK(k0,k1,N,dt,T,u0fun);

% Plot final time snapshot
figure;
plot(x, U(:,end),'r.' ,'LineWidth', 2);
xlabel('x'); ylabel('u(x,T)');
title(sprintf('Variable-k heat equation at T=%.3f', T));
grid on;

