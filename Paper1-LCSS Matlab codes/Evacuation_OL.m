clc
close all
clear all;

% Create a VideoWriter object
videoObj = VideoWriter('Evacuation_RK4_OL_alpha1_lambda0.1.mp4', 'MPEG-4');
videoObj.FrameRate = 10; % Adjust frame rate as needed
open(videoObj);

% Parameters
N = 500; % Number of individuals in the crowd
alpha1 = 1; % You can change the range and step size as needed
alpha2 = .2; % Repulsion strength between crowd and shooter
alpha3 = -1;
alpha4 = 5; % Strength of shooter's attraction towards crowd
alpha5 = 5; % Proportional control gain
alpha6 = .1;

lambda = 0.1;
p = 2.5; % Power-law decay of interactions
T = 4.195; % Total simulation time
dt = 0.05; % Time step
t = 0:dt:T;
num_steps = length(t);
rs = 0.4; % Define the distance threshold for killing a prey
rw = 0.5;
r0 = 4.7;
% Initial positions of individuals, shooter, and guided particle
% x0 = 2*rand(2, N); % Random initial positions for individuals
loadedData = load('x0.mat');
matrixStruct = loadedData.x0;
x0 = double(matrixStruct);
z0 = [3; 1]; % Initial position of shooter
y0 = [0.8; 0.8]; % Initial position of guided particle

% Room boundaries
room = [0, 5, 0, 5]; % [xmin, xmax, ymin, ymax]

target = .75*[room(2); room(4)]; % Top right corner of the room

killed_students = [];
num_runs = 1; % Number of runs for each "a" value to account for stochastic behavior
num_killed = 0;

% Time integration (RK4 method)

x = zeros(2, N, num_steps);
z = zeros(2, num_steps);
y = zeros(2, num_steps);
x(:, :, 1) = x0;
z(:, 1) = z0;
y(:, 1) = y0;

fig = figure();

for i = 1:num_steps
    
    k1 = motion(x(:,:,i), y(:,i), z(:,i),  N, room, target, p, lambda, rw, r0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6);
    k2 =  motion(x(:,:,i)+ 0.5*dt*k1(:,1:N), y(:,i)+ 0.5*dt*k1(:,N+1), z(:,i)+ 0.5*dt*k1(:,N+2), N, room, target, p, lambda, rw, r0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6);
    k3 =  motion(x(:,:,i)+ 0.5*dt*k2(:,1:N), y(:,i)+ 0.5*dt*k2(:,N+1), z(:,i)+ 0.5*dt*k2(:,N+2), N, room, target, p, lambda, rw, r0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6);
    k4 =  motion(x(:,:,i)+ 0.5*dt*k3(:,1:N), y(:,i)+ 0.5*dt*k3(:,N+1), z(:,i)+ 0.5*dt*k3(:,N+2), N, room, target, p, lambda, rw, r0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6);
    
    % Compute the next values using the weighted average of the intermediate values
    x(:,:,i+1) = x(:,:,i) + (dt/6)*(k1(:,1:N) + 2*k2(:,1:N) + 2*k3(:,1:N) + k4(:,1:N));
    y(:,i+1) = y(:,i) + (dt/6)*(k1(:,N+1) + 2*k2(:,N+1) + 2*k3(:,N+1) + k4(:,N+1));
    z(:,i+1) = z(:,i) + (dt/6)*(k1(:,N+2) + 2*k2(:,N+2) + 2*k3(:,N+2) + k4(:,N+2));
    
    %Terminate student is shot
    j = 1;
    
    while j <= N
        if  norm(x(:, j, i) - z(:, i)) < rs
            x = termination_procedure(x, j);
            N = N - 1; % Update the number of agents
            num_killed = num_killed + 1; % Increment the counter
        else
            j = j + 1;
        end
    end
    
    clf;
    
    % Draw the room walls
    fill([room(1)-.1, room(1),room(1),room(1)-.1 ],[room(3)-0.1, room(3)-0.1,room(4)+0.1, room(4)+0.1], 'k');  hold on;
    fill([room(2)+.1, room(2),room(2),room(2)+.1 ],[room(3)-0.1, room(3)-0.1,room(4)+0.1, room(4)+0.1], 'k');  % 'b' for blue color, 'FaceAlpha' for transparency
    fill([room(1)-.1, room(2)+.1,room(2)+0.1,room(1)-.1 ],[room(3), room(3),room(3)-0.1, room(3)-0.1], 'k');  % 'b' for blue color, 'FaceAlpha' for transparency
    fill([room(1)-.1, room(2)+.1,room(2)+0.1,room(1)-.1 ],[room(4), room(4),room(4)+0.1, room(4)+0.1], 'k');  % 'b' for blue color, 'FaceAlpha' for transparency
    
    % Shade the danger zone in red
    patch([room(1), 2.5, 2.5, room(1)], [room(3), room(3), 2.5, 2.5], [1, 0, 0], 'FaceAlpha', 0.2);
    
    % Shade the safe zone in green
    patch([2.5, room(2), room(2), 2.5], [2.5, 2.5, room(4), room(4)], [0, 1, 0], 'FaceAlpha', 0.2);
    
    alive_indices = find(~isnan(x(1, :, i)));
    scatter(x(1, alive_indices, i), x(2, alive_indices, i), 50, 'b', 'filled');
    hold on;
    plot(z(1, 1:i), z(2, 1:i),'-.r','linewidth',3);
    scatter(z(1, i+1), z(2, i+1),300,'rx', 'LineWidth', 200);
    
    plot(y(1, 1:i), y(2, 1:i), '-.g','linewidth',2);
    scatter(y(1, i+1), y(2, i+1),200, 'g+','linewidth',100);
    set(gca,'FontName','times','Fontsize',25)% yticklabels(G.Nodes(:,2)','Interpreter',"latex",'Fontname','times','Fontsize',15)
    
    
    % Display the information
    time_killed_text = sprintf('$t$: %.2f, $K$: %d', dt*i, num_killed);
    text( room(1) + 1.4, room(4) + 0.3,time_killed_text,'Interpreter',"latex",'Fontname','times','Fontsize',30, 'FontWeight', 'bold');
    xlim([room(1)-0.1, room(2)+0.1]);
    ylim([room(3)-.1, room(4)+0.1]);
    set(gcf, 'color', 'w');
    
    frame = getframe(gcf); % 'gcf' captures the current figure
    
    % Write the frame to the video file
    writeVideo(videoObj, frame);
    
    pause(0.001);
end


close(videoObj);% [min_killed_agents, min_idx] = min(killed_agents);


function dX = motion(x, y, z, N, room, target, p, lambda, rw, r0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6)
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
    dx(:, j) = dx(:, j) + alpha2 * (x(:, j) - z) / norm(x(:, j) - z)^2 ;
    dx(:, j) = dx(:, j) +alpha3* (x(:, j) - y);
    dx(:, j) = wall_constraint(dx(:, j), x(:, j), room, rw, alpha6);
end
%%%%
dz = zeros(size(z));
if z(1) <= 2.5 || z(2) <= 2.5
    for k = 1:N
        dz = dz + (x(:, k) - z) / norm(x(:, k) - z)^p;
    end
    dz = alpha4 / N * dz;
    
else
    % Terminate the simulation loop
    %     return
end

dz = wall_constraint(dz, z, room, rw, alpha6);
%%%%%

% Define the dynamics of the guiding particle

stop_threshold = 0.3; % Distance threshold for stopping

% Compute the weighted average position of the crowd
crowd_center = mean(x, 2);

% Compute the direction from the shooter to the crowd center
if z(1) <= 2.5 || z(2) <= 2.5
    away_from_shooter = crowd_center - z;
else
    away_from_shooter = 0;
end


desired_direction = lambda * (target - y) + (1 - lambda) * away_from_shooter;
% Check the distance to the target
distance_to_target = norm(target - y);



% Check the approaching speed to the target
%approaching_speed = dot(dy, (target - y)) / distance_to_target;

if distance_to_target < stop_threshold %&& approaching_speed < speed_threshold
    dy = [0; 0]; % Stop the guided particle
else
    % Calculate the velocity vector
    dy = alpha5 * desired_direction; % Proportional control law
    dy = wall_constraint(dy, y, room, rw, alpha6);
end
dX = [dx,dy,dz];
end


function x = termination_procedure(x, shot_index)
x(:, shot_index, :) = [];
end
function dp = wall_constraint(dp, p, room, rw, alpha6)

% Repulsive forces from walls
wall_fx = 0;
wall_fy = 0;

if p(1) < room(1) + rw
    wall_fx = alpha6 * (1/(p(1) - room(1)) - 1/rw);
elseif p(1) > room(2) - rw
    wall_fx = -alpha6 * (1/(room(2) - p(1)) - 1/rw);
end

if p(2) < room(3) + rw
    wall_fy = alpha6 * (1/(p(2) - room(3)) - 1/rw);
elseif p(2) > room(4) - rw
    wall_fy = -alpha6 * (1/(room(4) - p(2)) - 1/rw);
end

dp(1) = dp(1) + wall_fx;
dp(2) = dp(2) + wall_fy;

% Hard constraints for walls
dp(1) = max(min(dp(1), room(2) - p(1)), room(1) - p(1));
dp(2) = max(min(dp(2), room(4) - p(2)), room(3) - p(2));
end