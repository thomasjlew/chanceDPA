%	RISKDPA_MAIN - Solves a path planning problem with stochastics wk.
%                  and keeps the probability of hitting an obstacle below
%                  a threshold "delta".
% 
% Syntax:  riskDPA_main
%
% Inputs:
%    None
%
% Outputs:
%    Cost function for all states
%    Optimal policy
%    Final risk at timestep 0 using optimal policy
%
% Example: 
%    riskDPA_main
%
% Other m-files required:   (optional: generate_kernel.m)

% Subfunctions: none
% MAT-files required:   w_kernel_1_67_11x11.mat 
% Dataset required:     None - map automaticly generated

% References:
%   [1] M. Ono, M. Pavone, Y. Kuwata and J. Balaram, “Chance-constrained 
%       dynamic programming with application to risk-aware robotic space 
%       exploration,” Autonomous Robots, 2015.

% Author:   Thomas Lew
% email:    lewt@ethz.ch
% Website:  https://github.com/thomasjlew/
% November 2017; Last revision: 23-November-2017

%------------- BEGIN CODE --------------

clear all;
close all;
% clc;

%%  Initialize workspace: map, start and end, dynamics
disp('------------------------------------------------------------------');
disp('---------------------- Chance_DPA main script --------------------');
disp('------------------------------------------------------------------');
disp('- Initializing workspace with some map');
% Debugging booleans
B_PLOT = true;
B_SAVE_GIF = false;
B_PERFORMANCE = false;
B_PLOT_TRAJECTORY = true;
rng(1);         % Init. seed
digits(64);     % Increase precision
risk_str = [];  % For legend plotting

% Create binary map with obstacles
% Convention:    1 = no obstacle, 0 = obstacle present
MAP_WIDTH = 100;
MAP_HEIGHT = 100;
global map; map = ones(MAP_HEIGHT, MAP_WIDTH);  % No obstacles. 
obstacles = define_paper1_obstacles();          % Define obstacles
map = add_obstacles_to_map(map, obstacles);     % Add obstacles on the map

% Define start and end of the map
x0 = [5, 40];
xG = [95, 80];
[X,Y] = meshgrid(1:MAP_WIDTH,1:MAP_HEIGHT);
all_xks_in_map = [X(:),Y(:)];
    

% Plot the map + start/end
if B_PLOT
    figure(1); imshow(map, 'InitialMagnification', 300); 
    title('Map'); hold on; grid on; axis on;
    xlabel('x'); ylabel('y');
    plot(x0(1),x0(2), 'r+', 'LineWidth',3, 'MarkerSize', 10); 
    plot(xG(1),xG(2), 'g+', 'LineWidth',3, 'MarkerSize', 10);
    legend('Start', 'End', 'Location', 'NorthEast');
end

% Define dynamics and cost functions
disp('- Defining dynamics, variables and cost functions');
global N; N = 50;       % Nb steps for Dynamic Programming, see [1]
dk = 5;                 % Allowable input range, see [1]
u_space = allowable_inputs_dk(dk);
sigma = 1.67;           % Noise variance, see [1]
% w = round(sigma^2 * randn(2,40));   % For old approach, see "OLD CODE"
% [wY,wX] = meshgrid(-1:1,-1:1); w = [wX(:),wY(:)]';          
% load('w_kernel_1_67_5x5.mat');      % for weighted average by gaussian
% % load('w_kernel_1_67_11x11.mat'); % for weighted average by gaussian 
% KERNEL_WIDTH  = size(w_kernel,2);
% KERNEL_HEIGHT = size(w_kernel,1);

% Save cost-to-go gk values in a look-up table for faster processing
global ALPHA; ALPHA = 1e-5;
gk_values_for_uk = ALPHA * sqrt(u_space(1,:) .* u_space(1,:) + ...
                                u_space(2,:) .* u_space(2,:));

global fk; fk= @(xk, uk, wk) xk + uk + wk;         % Dynamics
global gN; gN = @(xN, xG) max(xN~=xG);             % Step cost
global gk; gk = @(xk, uk, alpha) alpha*norm(uk);   % Step cost at stage k
% Lk = @(k, xk, uk, lambda)         % Lagrangian step-wise (at script end)

% -------------------------------------------------------------------------
%%  Solve the problem using Chance-Constrained Dynamic Programming [1]
delta = 0.01;   % Probability constraint to be met using the optimal policy
e_tol = 1e-5;   % Error tolerance

% List of states to consider, might be changed for faster convergence
[Y,X] = meshgrid(1:MAP_HEIGHT, 1:MAP_WIDTH);
all_xks = [X(:),Y(:)];

% Gaussian filter to apply stochastics on each state (weighted average)
load('w_kernel_1_67_11x11.mat');
KERNEL_WIDTH = size(w_kernel,2);
KERNEL_HEIGHT = size(w_kernel,1);

% Save all optimal policies
best_mus = zeros([size(map),2, N+1]);

% -------------------------------------------------------------------------
% 1-5) Compute J_0^0(x0) (without constraint)
lambda = 0.0;
% lambda = 1e-5;

% -------------------------------------------------------------------------
% ----                  DYNAMIC PROGRAMMING ALGORITHM                  ----
% -------------------------------------------------------------------------
% K == N
% ------
% Save Ik values in a look-up table for faster processing
Ikxlambda_values_for_xk = zeros(size(map));
% Risk-to-go
r = zeros(size(map));
for xk=all_xks_in_map'
    Ikxlambda_values_for_xk(xk(2),xk(1)) = lambda * Ik(xk, map);
    r(xk(2),xk(1)) = Ik(xk, map);
    
    % Cost-to-go
    J(xk(2),xk(1)) = gN(xk', xG) + Ikxlambda_values_for_xk(xk(2),xk(1));
end
new_J = J;
best_mu = zeros([size(map),2]); % Best policy for each state
% best_mu_cur_state = xG.*ones(N+1,2);
best_mus(:,:,:,N+1) = best_mu;
% -------------------------------------------------------------------------


% -------------------
% K == [(N-1) ---> 0]
% -------------------
for k=N-1:-1:0
%     surf(r)
%     k
    
    % Apply stochastics as a convolution with noise kernel
    J_padded = padarray(J,[(KERNEL_HEIGHT-1)/2,(KERNEL_WIDTH-1)/2], ...
                        'replicate');
    J_w = conv2(J_padded, w_kernel, 'valid'); % removes padded edges
    
    % Obtain cost tensor with the value for each uk on the 3rd dim
    J_u_pad = repmat(padarray(J_w,[dk,dk],'replicate'),1,1,length(u_space));
    for uk_id = 1:length(u_space)
        uk = u_space(:,uk_id);
        J_u_pad(:,:,uk_id) = gk_values_for_uk(uk_id) + padarray(...
                padarray(J_w,[dk-uk(2),dk-uk(1)],'pre', 'replicate'), ...
                             [dk+uk(2),dk+uk(1)],'post','replicate');
                   
    end
    
    % Remove padded edges
    J_u = J_u_pad(dk+1:end-dk,dk+1:end-dk,:);
    
    % Get optimal policy and minimal cost value for each state
    [J,best_u_ids] = min(J_u,[],3); 
    for xk = all_xks'
        best_u_id = best_u_ids(xk(2), xk(1));
        best_mu(xk(2), xk(1), :) = u_space(:,best_u_id);
    end
     
    % Risk-to-go   (given best policy for this timestep)
    r_padded = padarray(r,[(KERNEL_HEIGHT-1)/2,(KERNEL_WIDTH-1)/2], ...
                        'replicate');
    r_w = conv2(r_padded, w_kernel, 'valid'); % also removes padded edges
    r_u_pad = repmat(padarray(r_w,[dk,dk],'replicate'),1,1,length(u_space));
    for uk_id = 1:length(u_space)
        uk = u_space(:,uk_id);
%         uk = [0;0]; %%%%%%%%
        r_u_pad(:,:,uk_id) = padarray(...
                padarray(r_w,[dk-uk(2),dk-uk(1)],'pre', 'replicate'), ...
                             [dk+uk(2),dk+uk(1)],'post','replicate');
                   
    end
    r_u = r_u_pad(dk+1:end-dk,dk+1:end-dk,:); % Remove padded edges
    for xk=all_xks_in_map'
        best_u_id = best_u_ids(xk(2), xk(1));
        r(xk(2),xk(1)) = Ik(xk, map) + r_u(xk(2),xk(1),best_u_id);
    end
    
    % Treat start separately (additionnal term if not k==0)
    if k ~= 0
        % Add lambda_k*Ik(xk)
        J = J + Ikxlambda_values_for_xk;
    end
    
    % Save next state if best policy & no stochastics
    best_mus(:,:,:,k+1) = best_mu;
    
    
    %% Plot results
    if B_SAVE_GIF
        set(gcf, 'Position', [800, 800, 800, 800])
        clf('reset') 
        figure(2)
        % subplot(2,1,1)
        % Cost J
        surf(J); hold on;
        plot_title = 'Cost J and optimal policy at timestep ' + string(k) + ...
                     ', lambda = ' + string(lambda);
        title(plot_title);
        xlabel('x'); ylabel('y'); zlabel('J')
        % Optimal policy
        [X,Y] = meshgrid(1:2:MAP_WIDTH,1:2:MAP_HEIGHT);
        mu_x_mat = best_mu(:,:,1); mu_y_mat = best_mu(:,:,2);
        quiver3(X,Y, J(1:2:end,1:2:end),best_mu(1:2:MAP_WIDTH,1:2:MAP_HEIGHT,1),...
                                 best_mu(1:2:MAP_WIDTH,1:2:MAP_HEIGHT,2), ...
                                 zeros(MAP_WIDTH/2,MAP_HEIGHT/2), 'r', 'LineWidth',1.5);
        % title('Optimal policy at at time step 0, lambda = 0');
        %xlabel('x'); ylabel('y');

        % End / Start
        scatter3(x0(1),x0(2), J(x0(2),x0(1)), 'c+', 'LineWidth',3, 'SizeData', 70); 
        scatter3(xG(1),xG(2), J(xG(2),xG(1)), 'g+', 'LineWidth',3, 'SizeData', 70);
        legend('Cost', 'Optimal policy', 'Start', 'End', 'Location', 'NorthWest');
        view([135 45])
        drawnow
        
        % Capture the plot as an image 
        filename = 'chanceDPA.gif'; 
        frame = getframe(2); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 

        % Write to the GIF File
        n = N-k
        if n == 1 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        end 
    end
end
% ----              END DYNAMIC PROGRAMMING ALGORITHM                  ----
% -------------------------------------------------------------------------

%% Get states reached with optimal policy assuming no stochastics
xk_s = zeros(N+1,2); xk_s(1,:) = x0;
for k = 1:N
    uk = best_mus(xk_s(k,2),xk_s(k,1),:,k);
    xk_s(k+1,:) = xk_s(k,:) + uk(:)';
end
if B_PLOT_TRAJECTORY
    figure(1);
    scatter(xk_s(:,1),xk_s(:,2),'LineWidth', 2);
    title('Optimal policies with collision risk probabilities');
    risk_str = [risk_str, ...
            'Opt. path, risk='+string(round(100*r(x0(2),x0(1)),3)) + '%'];
    legend('Start', 'End', risk_str(1),...
        'Location', 'NorthEastOutside');
    drawnow
end

% (NOT NECESSARY FOR ALGORITHM, ONLY FOR PLOTTING)
% Compute policy with interesting path for plotting
lambda=1e-4;
applyRiskDPA;
if B_PLOT_TRAJECTORY
    figure(1);
    scatter(xk_s(:,1),xk_s(:,2),'LineWidth', 2);
    risk_str = [risk_str, ...
            'Opt. path, risk='+string(round(100*r(x0(2),x0(1)),3)) + '%'];
    legend('Start', 'End', risk_str(1), risk_str(2), ...
        'Location', 'NorthEastOutside');
    drawnow
end
% (NOT NECESSARY FOR ALGORITHM, ONLY FOR PLOTTING)

%% Main DPA loop to obtain optimal policy 
% Skip step 2-8 (can be easily implemented ...)
% (contact me @ https://github.com/thomasjlew/chanceDPA/ if you need it)

lambda_L = 0;
lambda_U = 1;

% Big while loop: steps 10-21
iter = 0;
% while (lambda_U-lambda_L)*(r(x0(2),x0(1))-delta)>e_tol
% No Brent method yet, so solution is very close to chance constraint but
% not necessarily lower
while abs((lambda_U-lambda_L)*(r(x0(2),x0(1))-delta)) > e_tol
    iter = iter+1
    cur_eps = (lambda_U-lambda_L)*(r(x0(2),x0(1))-delta)
    
    lambda = (lambda_U-lambda_L)/2
    
    applyRiskDPA;
    
    if (r(x0(2),x0(1))-delta) == 0
        break
    elseif (r(x0(2),x0(1))-delta) < 0
        lambda_U = lambda;
    else
        lambda_L = lambda;
    end
end
disp('DONE!');
% End algorithm: the optimal policy meets the chance constraint.

B_PLOT_TRAJECTORY = true;
if B_PLOT_TRAJECTORY
    figure(1);
    scatter(xk_s(:,1),xk_s(:,2),'LineWidth', 2);
    risk_str = [risk_str, ...
            'Opt. path, risk='+string(round(100*r(x0(2),x0(1)),3)) + '%'];
    legend('Start', 'End', risk_str(1), risk_str(2), risk_str(3), ...
        'Location', 'NorthEastOutside');
    drawnow
end
        

% %% Plot results
% if B_PLOT
%     figure(2)
%     % subplot(2,1,1)
%     % Cost J
%     surf(J); hold on;
%     plot_title = 'Cost J and optimal policy at timestep ' + string(k) + ...
%                  ', lambda = ' + string(lambda);
%     title('plot_title');
%     xlabel('x'); ylabel('y'); zlabel('J')
%     % Optimal policy
%     [X,Y] = meshgrid(1:2:MAP_WIDTH,1:2:MAP_HEIGHT);
%     mu_x_mat = best_mu(:,:,1); mu_y_mat = best_mu(:,:,2);
%     quiver3(X,Y, J(1:2:end,1:2:end),best_mu(1:2:MAP_WIDTH,1:2:MAP_HEIGHT,1),...
%                              best_mu(1:2:MAP_WIDTH,1:2:MAP_HEIGHT,2), ...
%                              zeros(MAP_WIDTH/2,MAP_HEIGHT/2), 'r');
%     title('Optimal policy at at time step 0, lambda = 0');
%     xlabel('x'); ylabel('y');
% 
%     % End / Start
%     scatter3(x0(1),x0(2), J(x0(2),x0(1)), 'c+', 'LineWidth',2, 'SizeData', 70); 
%     scatter3(xG(1),xG(2), J(xG(2),xG(1)), 'g+', 'LineWidth',2, 'SizeData', 70);
%     legend('Cost', 'Optimal policy', 'Start', 'End');
% end

%% Debug functions commands
% Ik([20,65])
% Ik([10,10])

%% Additionnal functions
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Returns indicator random variable Ik (see [1])
% If xk is in the allowable states, return 0. Otherwise, returns 1
% 
% 23Nov: To minimize impact in efficiency, called once for each lambda to
%           initialize a look-up table.
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
% 23Nov: UNUSUED FUNCTION SINCE TOO SLOW IN MATLAB TO CALL THIS FUNCTION
% Returns value of step-wise Lagrangian
function Lk_val = Lk(k, xk, uk, lambda, map, ALPHA, N, xG, gk, gN, gk_values_for_uk, uk_id)
%     global ALPHA; %                       60% faster when removing global
%     global N;
%     global xG;
%     global gk;
%     global gN;
    
    if k==0
        Lk_val = gk_values_for_uk(uk_id);
%         Lk_val = ALPHA*(uk(1)*uk(1)+uk(2)+uk(2));%norm(uk);  %norm costly
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
% Function which returns similar obstacles as in [1]
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
    
    obstacle.pos = [78,80];   % x,y
    obstacle.size = [20, 20]; % width, height
    obstacles = [obstacles; obstacle];
end
%--------------------------------------------------------------------------


%% Old code
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% % Old DPA loop (slow because of for loops)
% Removed long loop, although it might be better for initialization...
% -------------------------------------------------------------------------
% for k=N-1:-1:N-4 
%     k
%     
%     if B_PERFORMANCE
%         profile on
%         tic;
%     end
%     
%     % Reduce state space since some states are useless to be updated
%     [Y,X] = meshgrid(max(1,xG(2)-(N-k)*dk-1):min(MAP_HEIGHT, xG(2)+(N-k)*dk+1), ...
%                      max(1,xG(1)-(N-k)*dk-1):min(MAP_WIDTH,  xG(1)+(N-k)*dk+1));
%     all_xks = [X(:),Y(:)];
%     
%     
%     
%     for xk = all_xks'
%         for uk_id = 1:length(u_space)
%             uk = u_space(:,uk_id);
%             % Compute reached state if no noise
%             xk_next_no_noise = xk+uk;
%             
%             % Check if out of bounds
%             if  max(xk_next_no_noise > [MAP_WIDTH-2; MAP_HEIGHT-2]) || ... % margin for wk
%                 max(xk_next_no_noise < [2; 2])
%                 J_u(uk_id) = MAX_J_INIT;
%                 continue
%             end
%             
%             % Advantage to reach goal
%             if xk_next_no_noise == xG'
%                 J_u = gN(xk_next_no_noise', xG) + gk_values_for_uk(uk_id);
%                 best_mu(xk(2), xk(1), :) = uk; 
%                 break
%             end
%             
%             % Compute Cost for this action, sampling for all noise values
%             w_nb_unitialized = 0;
%             J_w = w_kernel;
%             for wk = w
%                 % Apply dynamics
%     %             x_next = fk(xk, uk, wk); % Slooow
%                 x_next = xk + uk + wk; 
%                 
%                 % Special case: Cost function not defined here yet
%                 if J(x_next(2), x_next(1)) == MAX_J_INIT % Unexplored yet %%% POTENTIAL BUG!!!
%                     J_w(wk(2)+2,wk(1)+2) = 0;
%                     w_nb_unitialized = w_nb_unitialized + 1;
%                 else
%                     % ------------------------------------------------------
%                     % Compute value of Lk, faster (-60%) than calling function
%                     if k==0
%                         Lk_val = gk_values_for_uk(uk_id);
%                     elseif k == N
%                         Lk_val = gN(xk, xG) + lambda*Ik(xk, map);
%                     else
% %                         Lk_val = gk_values_for_uk(uk_id);% + lambda*Ik(xk, map);
%                         Lk_val = gk_values_for_uk(uk_id) + Ikxlambda_values_for_xk(xk(2),xk(1));
%                     end
% %                     J_w = J_w + Lk_val + J(x_next(2), x_next(1));
%                     % Matlab indices start at 1, HARDCODED FOR SPEED
%                     J_w(wk(2)+2,wk(1)+2) = J_w(wk(2)+2,wk(1)+2)*(Lk_val + J(x_next(2), x_next(1)));
%                     % ------------------------------------------------------
% 
%                     % ------------------------------------------------------
%                     % Call to Lk function (very very slow)
%     %                 J_w = J_w + Lk(k, xk, uk, lambda, map, ...
%     %                            ALPHA, N, xG, gk, gN, gk_values_for_uk, uk_id) + J(x_next(2), x_next(1));
%                     % ------------------------------------------------------
%                 end
%             end
%             
%             % Renormalize in case of unexplored states
%             if w_nb_unitialized ~= 9
%                 J_u(uk_id) = sum(sum(J_w))*(9/(9-w_nb_unitialized)); %%%% MAGIC NUMBER FOR KERNEL SIZE CHANGE THIS!!!
%             else
%                 J_u(uk_id) = MAX_J_INIT;
%             end
%         end
%         
%         % DPA: Select to minimum cost corresponding to best policy
%         [new_J(xk(2), xk(1)), best_u_id] = min(J_u);
%         best_mu(xk(2), xk(1), :) = u_space(:,best_u_id);
%     end
%     J = new_J; % don't forget this, or you'll obtain a Gauss-Seidel update
% 
%     if B_PERFORMANCE
%         toc
%     end
% end
% 
% disp('hey!');
% k
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% % Old DPA loop (slow because of for loops)
% load('w_kernel_1_67_5x5.mat'); % for weighted average by gaussian %%%%%
% J_u = [];
% for k=0:0
%     k
%     
%     if B_PERFORMANCE
%         profile on
%         tic;
%     end
%     
%     % Reduce state space since some states are useless to be updated
%     [Y,X] = meshgrid(max(1,xG(2)-(N-k)*dk-1):min(MAP_HEIGHT, xG(2)+(N-k)*dk+1), ...
%                      max(1,xG(1)-(N-k)*dk-1):min(MAP_WIDTH,  xG(1)+(N-k)*dk+1));
%     all_xks = [X(:),Y(:)];
%     
% %     all_xks = xG + allowable_inputs_dk((N-k)*dk+1)';
% %     all_xks = all_xks(min(ismember(all_xks,all_xks_in_map)')',:);
%     
%     % Remove states which were updated plenty times
% %     if k<=44
% %         updated_xks = xG + allowable_inputs_dk((N-(k+6))*dk)';
% %         all_xks = all_xks(max([ ~ismember(all_xks(:,1),updated_xks(:,1))';...
% %                                 ~ismember(all_xks(:,2),updated_xks(:,2))'])',:); 
%     
%     
%     
%     for xk = all_xks'
%         for uk_id = 1:length(u_space)
%             uk = u_space(:,uk_id);
%             % Compute reached state if no noise
%             xk_next_no_noise = xk+uk;
%             
%             % Check if out of bounds
%             if  max(xk_next_no_noise > [MAP_WIDTH-2; MAP_HEIGHT-2]) || ... % margin for wk
%                 max(xk_next_no_noise < [2; 2])
%                 J_u(uk_id) = MAX_J_INIT;
%                 continue
%             end
%             
%             % Advantage to reach goal
%             if xk_next_no_noise == xG'
%                 J_u = gN(xk_next_no_noise', xG) + gk_values_for_uk(uk_id);
%                 best_mu(xk(2), xk(1), :) = uk; 
%                 break
%             end
%             
%             % Compute Cost for this action, sampling for all noise values
%             w_nb_unitialized = 0;
%             J_w = w_kernel;
%             for wk = w
%                 % Apply dynamics
%     %             x_next = fk(xk, uk, wk); % Slooow
%                 x_next = xk + uk + wk; 
%                 
%                 % Special case: Cost function not defined here yet
%                 if J(x_next(2), x_next(1)) == MAX_J_INIT % Unexplored yet %%% POTENTIAL BUG!!!
%                     J_w(wk(2)+2,wk(1)+2) = 0;
%                     w_nb_unitialized = w_nb_unitialized + 1;
%                 else
%                     % ------------------------------------------------------
%                     % Compute value of Lk, faster (-60%) than calling function
%                     if k==0
%                         Lk_val = gk_values_for_uk(uk_id);
%                     elseif k == N
%                         Lk_val = gN(xk, xG) + lambda*Ik(xk, map);
%                     else
% %                         Lk_val = gk_values_for_uk(uk_id);% + lambda*Ik(xk, map);
%                         Lk_val = gk_values_for_uk(uk_id) + Ikxlambda_values_for_xk(xk(2),xk(1));
%                     end
% %                     J_w = J_w + Lk_val + J(x_next(2), x_next(1));
%                     % Matlab indices start at 1, HARDCODED FOR SPEED
%                     J_w(wk(2)+2,wk(1)+2) = J_w(wk(2)+2,wk(1)+2)*(Lk_val + J(x_next(2), x_next(1)));
%                     % ------------------------------------------------------
% 
%                     % ------------------------------------------------------
%                     % Call to Lk function (very very slow)
%     %                 J_w = J_w + Lk(k, xk, uk, lambda, map, ...
%     %                            ALPHA, N, xG, gk, gN, gk_values_for_uk, uk_id) + J(x_next(2), x_next(1));
%                     % ------------------------------------------------------
%                 end
%             end
%             
%             % Renormalize in case of unexplored states
%             if w_nb_unitialized ~= 9
%                 J_u(uk_id) = sum(sum(J_w))*(9/(9-w_nb_unitialized)); %%%% MAGIC NUMBER FOR KERNEL SIZE CHANGE THIS!!!
%             else
%                 J_u(uk_id) = MAX_J_INIT;
%             end
%         end
%         
%         % DPA: Select to minimum cost corresponding to best policy
%         [new_J(xk(2), xk(1)), best_u_id] = min(J_u);
%         best_mu(xk(2), xk(1), :) = u_space(:,best_u_id);
%     end
%     J = new_J; % don't forget this, or you'll obtain a Gauss-Seidel update
% 
%     if B_PERFORMANCE
%         toc
%     end
% end
% % END OLD DPA LOOP
%--------------------------------------------------------------------------

% END CODE
%--------------------------------------------------------------------------