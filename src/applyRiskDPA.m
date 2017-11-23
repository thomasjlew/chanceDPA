%	APPLYRISKDPA - DPA algorithm implementing the computations of the
%	optimal cost and expected risk of hitting an obstacle, as in [1].
% 
% Syntax:  main
%
% Inputs:
%    Parameters as initialized in riskDPA_main.
%
% Outputs:
%    Cost function for all states
%    Expected risk of hitting an obstacle following the optimal policy up
%       to timestep N
%
% Example: 
%    applyRiskDPA
%
% Other m-files required:   riskDPA_main and its dependencies

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
% November 2017; Last revision: 23-November-2017

%------------- BEGIN CODE --------------



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
%     k
    uk = best_mus(xk_s(k,2),xk_s(k,1),:,k);
    xk_s(k+1,:) = xk_s(k,:) + uk(:)';
end
% if B_PLOT_TRAJECTORY
%     scatter(xk_s(:,1),xk_s(:,2))
%     drawnow
% end


%--------------------------------------------------------------------------
% Returns indicator random variable Ik (see [1])
% If xk belongs to the allowable states, return 0. Otherwise, returns 1
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

% END CODE
%--------------------------------------------------------------------------