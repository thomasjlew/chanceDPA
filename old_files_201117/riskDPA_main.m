%	RISKDPA_MAIN - ...
% 
% Syntax:  main
%
% Inputs:
%    None
%
% Outputs:
%    Cost function for all states
%
% Example: 
%    main
%
% Other m-files required:   (optional: generate_kernel.m)

% Subfunctions: none
% MAT-files required:   w_kernel_... 
% Dataset required:     None

% References:
%   [1] M. Ono, M. Pavone, Y. Kuwata and J. Balaram, “Chance-constrained 
%       dynamic programming with application to risk-aware robotic space 
%       exploration,” Autonomous Robots, 2015.

% Author:   Thomas Lew
% email:    lewt@ethz.ch
% Website:  https://github.com/thomasjlew/
% November 2017; Last revision: 20-November-2017

%------------- BEGIN CODE --------------

clear all;
close all;
% clc;

%%  Initialize workspace: map, start and end, dynamics
disp('------------------------------------------------------------------');
disp('---------------------- Chance_DPA main script --------------------');
disp('------------------------------------------------------------------');
disp('Initializing workspace with some map');
% Debugging booleans
B_PLOT = true;
B_PERFORMANCE = false;
rng(1);

% Create binary map
% Convention:    1 = no obstacle, 0 = obstacle present
MAP_WIDTH = 100;
MAP_HEIGHT = 100;
global map;
map = ones(MAP_HEIGHT, MAP_WIDTH);     % no obstacles. 

% Define and add obstacles on the map
obstacles = define_paper1_obstacles();
map = add_obstacles_to_map(map, obstacles);

% Define start and end on the map
x0 = [5, 40];
global xG;
xG = [95, 80];
[X,Y] = meshgrid(1:MAP_WIDTH,1:MAP_HEIGHT);
all_xks_in_map = [X(:),Y(:)];
    

% Plot the map + start/end
if B_PLOT
    figure(1); imshow(map, 'InitialMagnification', 600); 
    title('Map'); hold on; grid on; axis on;
    xlabel('x'); ylabel('y');
    scatter(x0(1),x0(2), 'r+', 'LineWidth',2, 'SizeData', 30); 
    scatter(xG(1),xG(2), 'g+', 'LineWidth',2, 'SizeData', 30);
end

% Define dynamics and cost functions
global N; N = 50;      % Nb steps for Dynamic Programming, see [1]
dk = 5;             % Allowable input range, see [1]
sigma = 1.67;          % Noise variance, see [1]
w = round(sigma^2 * randn(2,40));
w = [0,0];                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% removed random thing !!!
w = [-1,-1];  
[wY,wX] = meshgrid(-1:1,-1:1); w = [wX(:),wY(:)]';          
load('w_kernel_1_67_5x5.mat'); % for weighted average by gaussian
% load('w_kernel_1_67_11x11.mat'); % for weighted average by gaussian %%%%%
% CHECK DIMENSIONS IF YOU CHANGE SOMETHING HERE !!!

global ALPHA; ALPHA = 1e-5;
u_space = allowable_inputs_dk(dk);
% Save gk values in a look-up table for faster processing
gk_values_for_uk = ALPHA * sqrt(u_space(1,:) .* u_space(1,:) + ...
                                u_space(2,:) .* u_space(2,:));

global fk; fk= @(xk, uk, wk) xk + uk + wk;         % Dynamics
global gN; gN = @(xN, xG) max(xN~=xG);             % Step cost
global gk; gk = @(xk, uk, alpha) alpha*norm(uk);   % Step cost at stage k
% Lk = @(k, xk, uk, lambda)         % Lagrangian step-wise (at script end)

% -------------------------------------------------------------------------
%%  Solve the problem using Chance-Constrained Dynamic Programming [1]
e_tol = 0.01; % Try different values ---< -< -<-- -<- -<--< ---------------------------------!!!
delta = 0.01;

% -------------------------------------------------------------------------
% 1-5) Compute J_0^0(x0) (without constraint)
lambda = 0;
% Save Ik values in a look-up table for faster processing
Ikxlambda_values_for_xk = zeros(size(map));
for xk=all_xks_in_map'
    Ikxlambda_values_for_xk(xk(2),xk(1)) = lambda * Ik(xk, map);
    J(xk(2),xk(1)) = gN(xk', xG) + Ikxlambda_values_for_xk(xk(2),xk(1));
end
% Initialize costs
MAX_J_INIT = 100;
% J = MAX_J_INIT*ones(size(map));
% J(xG(2), xG(1)) = Lk(N, xG, [], lambda, map, ALPHA, N, xG, gk, gN, gk_values_for_uk);

new_J = J;
best_mu = zeros([size(map),2]); % Best policy for each state

tic
profile on
for k=N-1:-1:0
    k
    
    if B_PERFORMANCE
        profile on
        tic;
    end
    
    % Reduce state space since some states are useless to be updated
    [Y,X] = meshgrid(max(1,xG(2)-(N-k)*dk-1):min(MAP_HEIGHT, xG(2)+(N-k)*dk+1), ...
                     max(1,xG(1)-(N-k)*dk-1):min(MAP_WIDTH,  xG(1)+(N-k)*dk+1));
    all_xks = [X(:),Y(:)];
    
%     all_xks = xG + allowable_inputs_dk((N-k)*dk+1)';
%     all_xks = all_xks(min(ismember(all_xks,all_xks_in_map)')',:);
    
    % Remove states which were updated plenty times
%     if k<=44
%         updated_xks = xG + allowable_inputs_dk((N-(k+6))*dk)';
%         all_xks = all_xks(max([ ~ismember(all_xks(:,1),updated_xks(:,1))';...
%                                 ~ismember(all_xks(:,2),updated_xks(:,2))'])',:); 
%     end
    
    if size(all_xks,1) == 0
        break % Everything is updated plenty of times!!! 
    end
    toc
    
    for xk = all_xks'
        for uk_id = 1:length(u_space)
            uk = u_space(:,uk_id);
            
            % Check if out of bounds
            xk_next_no_noise = xk+uk;
            if  max(xk_next_no_noise > [MAP_WIDTH-2; MAP_HEIGHT-2]) || ... % margin for wk
                max(xk_next_no_noise < [2; 2])
                J_u(uk_id) = MAX_J_INIT;
                continue
            end
            
            % Advantage to reach goal
            if xk_next_no_noise == xG'
                J_u = gN(xk_next_no_noise', xG) + gk_values_for_uk(uk_id);
                best_mu(xk(2), xk(1), :) = uk; 
                break
            end
            
            % Compute Cost for this action, sampling for all noise values
            w_nb_unitialized = 0;
            J_w = w_kernel;
            for wk = w
                % Apply dynamics
    %             x_next = fk(xk, uk, wk); % Slooow
                x_next = xk + uk + wk; 
                
                % Special case: Cost function not defined here yet
                if J(x_next(2), x_next(1)) == MAX_J_INIT % Unexplored yet %%% POTENTIAL BUG!!!
                    J_w(wk(2)+2,wk(1)+2) = 0;
                    w_nb_unitialized = w_nb_unitialized + 1;
                else
                    % ------------------------------------------------------
                    % Compute value of Lk, faster (-60%) than calling function
                    if k==0
                        Lk_val = gk_values_for_uk(uk_id);
                    elseif k == N
                        Lk_val = gN(xk, xG) + lambda*Ik(xk, map);
                    else
%                         Lk_val = gk_values_for_uk(uk_id);% + lambda*Ik(xk, map);
                        Lk_val = gk_values_for_uk(uk_id) + Ikxlambda_values_for_xk(xk(2),xk(1));
                    end
%                     J_w = J_w + Lk_val + J(x_next(2), x_next(1));
                    % Matlab indices start at 1, HARDCODED FOR SPEED
                    J_w(wk(2)+2,wk(1)+2) = J_w(wk(2)+2,wk(1)+2)*(Lk_val + J(x_next(2), x_next(1)));
                    % ------------------------------------------------------

                    % ------------------------------------------------------
                    % Call to Lk function (very very slow)
    %                 J_w = J_w + Lk(k, xk, uk, lambda, map, ...
    %                            ALPHA, N, xG, gk, gN, gk_values_for_uk, uk_id) + J(x_next(2), x_next(1));
                    % ------------------------------------------------------
                end
            end
            
            % Renormalize in case of unexplored states
            if w_nb_unitialized ~= 9
                J_u(uk_id) = sum(sum(J_w))*(9/(9-w_nb_unitialized)); %%%% MAGIC NUMBER FOR KERNEL SIZE CHANGE THIS!!!
            else
                J_u(uk_id) = MAX_J_INIT;
            end
        end
        
        % DPA: Select to minimum cost corresponding to best policy
        [new_J(xk(2), xk(1)), best_u_id] = min(J_u);
        best_mu(xk(2), xk(1), :) = u_space(:,best_u_id);
    end
    J = new_J;

    if B_PERFORMANCE
        toc
    end
end
% for k=40:N
%     for xk = reachable_states(
%     
%     J(
% end
% -------------------------------------------------------------------------
toc

% Plot result
figure
surf(J)
title('Cost J at time step 0, lambda = 0');
xlabel('x'); ylabel('y'); zlabel('J')

%% Debug functions commands
% Ik([20,65])
% Ik([10,10])

%% Additionnal functions
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Returns indicator random variable Ik (see [1])
% If xk is in the allowable states, return 0. Otherwise, returns 1
function Ik_val = Ik(xk, map)
%    global map;
   
   if map(xk(2),xk(1)) ==  0
       Ik_val = 1;
       return
   else
       Ik_val = 0;
       return
   end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Returns value of step-wise Lagrangian
function Lk_val = Lk(k, xk, uk, lambda, map, ALPHA, N, xG, gk, gN, gk_values_for_uk, uk_id)
%     global ALPHA; %                       60% faster when removing global
%     global N;
%     global xG;
%     global gk;
%     global gN;
    
    if k==0
        Lk_val = gk_values_for_uk(uk_id);
%         Lk_val = ALPHA*(uk(1)*uk(1)+uk(2)+uk(2));%norm(uk);             %%% REMOVEDE THE NORM WHICH TAKES WAY TOO MUCH TIME
%         Lk_val = gk(xk, uk, ALPHA);
        return;
    elseif k == N
        Lk_val = gN(xk, xG) + lambda*Ik(xk, map);
        return;
    else
        Lk_val = gk_values_for_uk(uk_id);
%         Lk_val = ALPHA*(uk(1)*uk(1)+uk(2)+uk(2));%norm(uk);
%         Lk_val = gk(xk, uk, ALPHA) + lambda*Ik(xk, map);
        return
    end
    
        
   
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Returns allowable inputs [2xN] which norm is smaller than "dist"
function u = allowable_inputs_dk(dist)
    u = [];
    
    [X,Y] = meshgrid(-dist:dist,-dist:dist); 
    x = X(:); y = Y(:);
    
    for i = 1:length(x)
        if norm([x(i); y(i)]) <= dist
            u = [u, [x(i); y(i)]];
        end
    end
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Adds an obstacle
function new_map = add_obstacle_to_map(new_map, obstacle)
    new_map(obstacle.pos(2)-round(obstacle.size(2)/2) : ...
            obstacle.pos(2)+round(obstacle.size(2)/2) , ...
            obstacle.pos(1)-round(obstacle.size(1)/2) : ...
            obstacle.pos(1)+round(obstacle.size(1)/2)) = 0;%...
%                                                 ones(obstacle.size(2), ...
%                                                      obstacle.size(1));
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Function to add multiple obstacles
function new_map = add_obstacles_to_map(new_map, obstacles)
    nb_obstacles = size(obstacles);
    
    for i=1:nb_obstacles        
        new_map = add_obstacle_to_map(new_map, obstacles(i));
    end
end
%--------------------------------------------------------------------------
        
%--------------------------------------------------------------------------
% Function which returns similara obstacles as in [1]
function obstacles = define_paper1_obstacles()
    obstacles = [];
    
    obstacle.pos = [40,18];   % x,y
    obstacle.size = [20, 20]; % width, height
    obstacles = [obstacles; obstacle];
    
    obstacle.pos = [45,45];   % x,y
    obstacle.size = [20, 20]; % width, height
    obstacles = [obstacles; obstacle];
    
    obstacle.pos = [20,65];   % x,y
    obstacle.size = [20, 20]; % width, height
    obstacles = [obstacles; obstacle];
    
    obstacle.pos = [48,70];   % x,y
    obstacle.size = [20, 20]; % width, height
    obstacles = [obstacles; obstacle];
    
    obstacle.pos = [78,85];   % x,y
    obstacle.size = [20, 20]; % width, height
    obstacles = [obstacles; obstacle];
end
%--------------------------------------------------------------------------


