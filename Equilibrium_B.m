clc
close all
clear all;

% Parameters
N = 500; % Number of individuals in the crowd
alpha1 = 1; % Repulsion/attraction strength between individuals
alpha2 = 0.2; % Repulsion strength between crowd and shooter
alpha3 = 0; % Attraction strength toward guiding agent
alpha4 = 5; % Strength of shooter's attraction towards crowd
lambda = 0.1; % Guiding agent's direction parameter
p = 2; % Power-law decay of interactions
T = 30; % Total simulation time
dt = 0.5; % Time step
t = 0:dt:T;
num_steps = length(t);
r0 = 4.7; % Decay effect parameter


R1 = sqrt(alpha2/(alpha1-alpha3)); % Radius of inner circle
R2 = sqrt((1+alpha2)/(alpha1-alpha3)); % Radius of outer circle

% Initial positions of individuals, shooter, and guiding agent
x0 = rand(2, N); % Random initial positions for individuals
y0 = [0.5; 0.5]; % Initial position of guiding agent
z0 = [0.5; 0.5]; % Initial position of shooter

theta = linspace(0, 2*pi, 100); % You can change the number of points (100 in this case)
xr1 = R1 * cos(theta);
yr1 = R1 * sin(theta);
xr2 = R2 * cos(theta);
yr2 = R2 * sin(theta);
x = zeros(2, N, num_steps);
z = zeros(2, num_steps);
y = zeros(2, num_steps);
x(:, :, 1) = x0;
z(:, 1) = z0;
y(:, 1) = y0;

for i = 1:num_steps
    
    % Compute intermediate values
    k1 = motion(x(:,:,i), y(:,i), z(:,i), N, p, r0, alpha1, alpha2,alpha3,alpha4);
    k2 =  motion(x(:,:,i)+ 0.5*dt*k1(:,1:N), y(:,i)+ 0.5*dt*k1(:,N+1), z(:,i)+ 0.5*dt*k1(:,N+2), N, p, r0, alpha1, alpha2,alpha3,alpha4);
    k3 =  motion(x(:,:,i)+ 0.5*dt*k2(:,1:N), y(:,i)+ 0.5*dt*k2(:,N+1), z(:,i)+ 0.5*dt*k2(:,N+2), N, p, r0, alpha1, alpha2,alpha3,alpha4);
    k4 =  motion(x(:,:,i)+ 0.5*dt*k3(:,1:N), y(:,i)+ 0.5*dt*k3(:,N+1), z(:,i)+ 0.5*dt*k3(:,N+2), N, p, r0, alpha1, alpha2,alpha3,alpha4);
    
    % Compute the next values using the weighted average of the intermediate values
    x(:,:,i+1) = x(:,:,i) + (dt/6)*(k1(:,1:N) + 2*k2(:,1:N) + 2*k3(:,1:N) + k4(:,1:N));
    y(:,i+1) = y(:,i) + (dt/6)*(k1(:,N+1) + 2*k2(:,N+1) + 2*k3(:,N+1) + k4(:,N+1));
    z(:,i+1) = z(:,i) + (dt/6)*(k1(:,N+2) + 2*k2(:,N+2) + 2*k3(:,N+2) + k4(:,N+2));
end


fig = figure();
scatter(x(1, :, end), x(2, :, end), 50, 'b', 'filled');
hold on;
centerx = (alpha1*(sum(x(1,:,num_steps))/N) - alpha3*y(1,end))/(alpha1-alpha3);
centery = (alpha1*(sum(x(2,:,num_steps))/N) - alpha3*y(2,end))/(alpha1-alpha3);
plot(xr1+z(1, end), yr1+z(2, end), 'k', 'LineWidth', 5);
plot(xr2+centerx, yr2+centery, 'k', 'LineWidth', 5);
plot(z(1, end), z(2, end), 'rx', 'MarkerSize', 30, 'LineWidth', 30);
set(gca,'color', 'w','FontName','times','Fontsize',17);
axis equal
axis off
saveas(fig,'caseB.eps')
saveas(fig,'caseB.fig')
saveas(fig,'caseB.jpg')

function dX = motion(x, y, z, N, p, r0, alpha1, alpha2,alpha3,alpha4)
dx = zeros(size(x));
for j = 1:N
    for k = 1:N
        if k ~= j
            r = norm(x(:, j) - x(:, k));
            if r <= r0
                h = 1;
            else
                h = exp(-(r-r0));
            end
            dx(:, j) = dx(:, j) + (x(:, j) - x(:, k)) / norm(x(:, j) - x(:, k))^2 - alpha1 * (x(:, j) - x(:, k))*h ; % Corrected force calculation
        end
    end
    dx(:, j) = (1/N)*dx(:, j);
    dx(:, j) = dx(:, j) + alpha2 * (x(:, j) - z) / norm(x(:, j) - z)^2 +alpha3* (x(:, j) - y);
end

%%%%
dz = zeros(size(z));
if z(1) <= 2.5 || z(2) <= 2.5
    for k = 1:N
        dz = dz + (x(:, k) - z) / norm(x(:, k) - z)^p;
    end
    dz = alpha4 / N * dz;
    
end
    
dy = [0; 0];

dX = [dx,dy,dz];
end



