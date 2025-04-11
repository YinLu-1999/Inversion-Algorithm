function [x_estimates,y_estimates,z_estimates,Mx_estimates,My_estimates,Mz_estimates] = text_inverse(x_points,y_points,d,x1,y1,Bz,dBxdx,dBxdy,dBxdz,dBydy,dBydz)
mu0 = 4*pi*1e-7;
% 初始化结果数组
    num_x = length(x_points);
    num_y = length(y_points);
    M_xyz_results = zeros(3, []);
    r_dipole_results = zeros(3, []);
    errors = [];
    
%计算每个区域每个点的一阶梯度矩阵
for ix = 1:num_x
    for iy = 1:num_y
        r_obs = [x_points(ix); y_points(iy); d];
        % 当前观测点位置
        % 找到最接近(0.5, 0.5)的网格点
       
        target_x = x_points(ix);
        target_y = y_points(iy);

% 找到最接近的x和y索引
[~, x_idx] = min(abs(x1 - target_x));
[~, y_idx] = min(abs(y1 - target_y));

        
        B_or = Bz(x_idx, y_idx);
        %当前张量矩阵
        dBxdx_current = dBxdx(x_idx, y_idx);
        dBxdy_current = dBxdy(x_idx, y_idx);
        dBxdz_current = dBxdz(x_idx, y_idx);
        dBydy_current = dBydy(x_idx, y_idx);
        dBydz_current = dBydz(x_idx, y_idx);
       
        % 构建张量矩阵
        D1 = [dBxdx_current,dBxdy_current,dBxdz_current; dBxdy_current,dBydy_current,dBydz_current; dBxdz_current,dBydz_current,-dBxdx_current-dBydy_current];
        
        % 计算特征值和特征向量
        [V,eigenvalues] = eig(D1);
        eigenvalues = diag(eigenvalues);
        
        % 保存原始的特征值和特征向量
        original_lambda = eigenvalues;
        original_V = V;
        
        % 生成所有可能的排列
        combinations = perms(eigenvalues);
        
        % 创建条件检查函数(n1不等于0)
        condition1 = @(x)( 2*x(2) + x(1) )/(x(2)-x(1))>= 0;
        condition2 = @(x) (2*x(2) + x(1))/(x(2)-x(1)) <= 1;
        
        % 找到满足所有条件的排列
        valid_combinations = [];
        valid_indices = {};
        
        for i = 1:size(combinations, 1)
            current = combinations(i,:);
            if condition1(current) && condition2(current)
                indices = zeros(1, length(current));
                for j = 1:length(current)
                    possible_indices = find(abs(original_lambda - current(j)) < 1e-10);
                    used_indices = indices(1:j-1);
                    available_indices = setdiff(possible_indices, used_indices);
                    if ~isempty(available_indices)
                        indices(j) = available_indices(1);
                    end
                end
                valid_combinations = [valid_combinations; current];
                valid_indices{end+1} = indices;
            end
        end
        
        % 处理满足条件的排列
        if ~isempty(valid_combinations)
            scores = valid_combinations(:,1) + 2*valid_combinations(:,2);
            [~, idx] = sort(scores, 'descend');
            final_lambda = valid_combinations(idx(1),:);
            final_indices = valid_indices{idx(1)};
            final_V = V(:,final_indices);
            
            % 计算参数
            lambda_1 = final_lambda(1);
            lambda_2 = final_lambda(2);
            lambda_3 = final_lambda(3);
            
            sin_theta_squared = (2*lambda_2 + lambda_1) / (lambda_2 - lambda_1);
            cos_theta_squared = 1 - sin_theta_squared;
            
            % 四种可能的组合
            sin_thetas = [sqrt(sin_theta_squared), -sqrt(sin_theta_squared), 
                        sqrt(sin_theta_squared), -sqrt(sin_theta_squared)];
            cos_thetas = [sqrt(cos_theta_squared), sqrt(cos_theta_squared), 
                        -sqrt(cos_theta_squared), -sqrt(cos_theta_squared)];
            
            % 存储每种组合的结果
            combo_errors = zeros(1, 4);
            combo_M_xyz = zeros(3, 4);
            combo_r_xyz = zeros(3, 4);
            
            % 遍历所有组合
            for combo = 1:4
                sin_theta = sin_thetas(combo);
                cos_theta = cos_thetas(combo);
                
                n = [cos_theta; sin_theta; 0];
                n1 = n(1);
                n2 = n(2);
                
                % 计算磁矩分量
                M1 = lambda_1 * n1 * (5 * n2^2 - 1) / (1 + n1^2);
                M2 = lambda_1 * n2 * (1 - 5 * n1^2) / (1 + n1^2);
                M3 = 0;
                M_calc = [M1; M2; M3];
                
                % 坐标转换
                Q = inv(final_V);
                alpha = Q(1,:);
                belta = Q(2,:);
                gamma = Q(3,:);
                
                % 归一化
                f_norm = sqrt(alpha(1)^2 + belta(1)^2 + gamma(1)^2);
                s_norm = sqrt(alpha(2)^2 + belta(2)^2 + gamma(2)^2);
                t_norm = sqrt(alpha(3)^2 + belta(3)^2 + gamma(3)^2);
                
                A = [alpha(1)/f_norm, alpha(2)/s_norm, alpha(3)/t_norm;
                     belta(1)/f_norm, belta(2)/s_norm, belta(3)/t_norm;
                     gamma(1)/f_norm, gamma(2)/s_norm, gamma(3)/t_norm];
                
                % 转换到xyz坐标系
                n_xyz = A\n;
                M_xyz = A\M_calc;
                
                % 计算 r
                r = 3 * B_or / (((3 * M_xyz(1) * n_xyz(1) + 3 * M_xyz(2) * n_xyz(2) + 3 * M_xyz(3) * n_xyz(3)) * n_xyz(3)) - M_xyz(3));
                
                % 计算位置向量
                r_xyz = r * n_xyz;
                M_xyz = (4*pi*r^4)*M_xyz/(3*mu0);
                
                % 计算磁场并比较误差
                r_mag = norm(r_xyz);
                B_calculated = (mu0/(4*pi)) * (3*(dot(M_xyz,r_xyz)*r_xyz)/(r_mag^5) - M_xyz/(r_mag^3));
                combo_errors(combo) = norm(B_calculated(3) - B_or);
                combo_M_xyz(:,combo) = M_xyz;
                combo_r_xyz(:,combo) = r_xyz;
            end
            
            % 选择最佳组合
            [min_error, best_combo] = min(combo_errors);
            
            % 存储结果
            M_xyz_results(:,ix,iy) = combo_M_xyz(:,best_combo);
            r_dipole_results(:,ix,iy) = r_obs - combo_r_xyz(:,best_combo);
            errors(ix,iy) = min_error;
        else
            % 如果没有找到有效解，标记为无效结果
            M_xyz_results(:,ix,iy) = NaN(3,1);
            r_dipole_results(:,ix,iy) = NaN(3,1);
            errors(ix,iy) = Inf;
        end
    end
end

%% 结果可视化
% 找到误差最小的点
[min_error_val, min_idx] = min(errors(:));
[min_ix, min_iy] = ind2sub(size(errors), min_idx);

% 位置

% 将所有估计位置转换为可绘制的格式
x_estimates = reshape(r_dipole_results(1,:,:), [], 1);
y_estimates = reshape(r_dipole_results(2,:,:), [], 1);
z_estimates = reshape(r_dipole_results(3,:,:), [], 1);

%磁矩
% 将所有估计位置转换为可绘制的格式
Mx_estimates = reshape(M_xyz_results(1,:,:), [], 1);
My_estimates = reshape(M_xyz_results(2,:,:), [], 1);
Mz_estimates = reshape(M_xyz_results(3,:,:), [], 1);
end