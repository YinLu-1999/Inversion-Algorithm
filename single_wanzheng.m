clc;clear;
close all;
%% 正向计算
% 参数值
r_dipole = [8; 5;2]; 
M = [16;32;66];
r_obs = [10;9;7]; 
mu0 =4*pi*1e-7;
% 测距
r_for = r_obs - r_dipole; % 从偶极子到观测点的向量
r_mag = norm(r_for);   % 距离大小
% B = (mu0 / (4 * pi * r_mag^3)) * (3 * (M' * r) * r - M * r_mag^2);
B = (mu0/(4*pi)) * (3*(dot(M,r_for)*r_for)/(r_mag^5) - M/(r_mag^3));
% 导数
x = r_for(1);
y = r_for(2);
z = r_for(3);
r_mag = norm(r_for);
constant = (3*mu0) / (4 * pi * r_mag^7);
dBxdx = constant * ((3*r_mag^2 - 5*x^2)* x * M(1) + (r_mag^2-5*x^2) * y * M(2)+(r_mag^2-5*x^2)*z*M(3));
dBxdy = constant * ((r_mag^2 - 5*x^2) * y* M(1) + (r_mag^2-5*y^2) * x * M(2)+(-5)*x*y*z*M(3));
dBxdz = constant * ((r_mag^2 - 5*x^2) * z* M(1) + (-5)*x* y*z * M(2)+(r_mag^2-5*z^2)*x*M(3));
dBydz = constant * ((-5)*x*y*z * M(1) + (r_mag^2-5*y^2) * z * M(2)+(r_mag^2-5*z^2)*y*M(3));
dBydy = constant * ((r_mag^2 - 5*y^2)*x * M(1) + (3*r_mag^2-5*y^2) * y * M(2)+(r_mag^2-5*y^2)*z*M(3));
%% 构建张量矩阵

D1 = [dBxdx,dBxdy,dBxdz; dBxdy, dBydy,dBydz;dBxdz,dBydz,-dBxdx-dBydy];
%% 计算特征值和特征向量
[V,eigenvalues] = eig(D1);
eigenvalues=diag(eigenvalues);



% 保存原始的特征值和特征向量
original_lambda = eigenvalues;
original_V = V;
% 生成所有可能的排列
combinations = perms(eigenvalues);

% % 创建条件检查函数(n1=0)
% condition1 = @(x)( x(2) - x(1) )/(x(1)+2*x(2))>= 0;  % 第一个条件
% condition2 = @(x) (x(2) - x(1))/(x(1)+2*x(2)) <= 1;    % 第二个条件

% 创建条件检查函数(n1不等于0)
condition1 = @(x)( 2*x(2) + x(1) )/(x(2)-x(1))>= 0;  % 第一个条件
condition2 = @(x) (2*x(2) + x(1))/(x(2)-x(1)) <= 1;    % 第二个条件

% 找到满足所有条件的排列
valid_combinations = [];
valid_indices = {};  % 用于存储每个有效排列对应的索引

for i = 1:size(combinations, 1)
    current = combinations(i,:);
    if condition1(current) && condition2(current)
        % 找到当前排列中每个值在原始特征值中的位置
        indices = zeros(1, length(current));
        for j = 1:length(current)
            % 找到当前值在原始特征值中的位置
            % 注意处理可能的重复特征值
            possible_indices = find(abs(original_lambda - current(j)) < 1e-10);
            % 选择还未使用的索引
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

% 如果找到了满足条件的排列
if ~isempty(valid_combinations)
    % 按照条件1的值大小排序
    scores = valid_combinations(:,1) + 2*valid_combinations(:,2);
    [~, idx] = sort(scores, 'descend');
    
    % 获取得分最高的排列及其对应的索引
    final_lambda = valid_combinations(idx(1),:);
    final_indices = valid_indices{idx(1)};
    
    % 重排特征向量
    final_V = V(:,final_indices);
    
    % 验证特征值和特征向量的对应关系
    for i = 1:length(final_lambda)
        residual = norm(D1*final_V(:,i) - final_lambda(i)*final_V(:,i));
        if residual > 1e-10
            warning(['特征值和特征向量的对应关系验证失败，第 ' num2str(i) ' 对的残差为 ' num2str(residual)]);
        end
    end
 % 输出结果
    disp('满足条件的特征值排序为：');
    disp(final_lambda);
    disp('对应的特征向量为：');
    disp(final_V);
else
    disp('没有找到满足所有条件的排列');
end


% 按原矩阵顺序输出特征值
lambda_1 = final_lambda(1);%算的时候用
lambda_2 = final_lambda(2);%
lambda_3 = final_lambda(3);%方向为0是最大的  
% ... [前面的代码保持不变直到计算 sin_theta_squared] ...

% 计算四种可能的组合
sin_theta_squared = (2*lambda_2 + lambda_1) / (lambda_2 - lambda_1);
cos_theta_squared = 1 - sin_theta_squared;

% 四种可能的组合
sin_thetas = [sqrt(sin_theta_squared), -sqrt(sin_theta_squared), 
              sqrt(sin_theta_squared), -sqrt(sin_theta_squared)];
cos_thetas = [sqrt(cos_theta_squared), sqrt(cos_theta_squared), 
              -sqrt(cos_theta_squared), -sqrt(cos_theta_squared)];

% 存储每种组合的结果和误差
num_combinations = 4;
errors = zeros(1, num_combinations);
all_M_xyz = zeros(3, num_combinations);
all_r_xyz = zeros(3, num_combinations);

% 遍历所有组合
for combo = 1:num_combinations
    % 当前组合的 sin_theta 和 cos_theta
    sin_theta = sin_thetas(combo);
    cos_theta = cos_thetas(combo);
    
    % 计算 n 向量
    n = [cos_theta; sin_theta; 0]
    n1 = n(1);
    n2 = n(2);
    n3 = n(3);
    
    % 计算磁矩分量
    M1 = lambda_1 * n1 * (5 * n2^2 - 1) / (1 + n1^2);
    M2 = lambda_1 * n2 * (1 - 5 * n1^2) / (1 + n1^2);
    M3 = 0;
    M = [M1; M2; M3];
    
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
    M_xyz = A\M
    
    % 计算 r
    Mx = M_xyz(1);
    My = M_xyz(2);
    Mz = M_xyz(3);
    nx = n_xyz(1);
    ny = n_xyz(2);
    nz = n_xyz(3);

    Bx=B(1);
    By=B(2);
    Bz=B(3);
    
    % 使用三个分量都计算r，取平均值增加稳定性
    r = 3 * Bz / (((3 * Mx * nx + 3 * My * ny + 3 * Mz * nz) * nz) - Mz);
    % r = 3 * Bx / (((3 * Mx * nx + 3 * My * ny + 3 * Mz * nz) * nx) - Mx);
    % r = 3 * By / (((3 * Mx * nx + 3 * My * ny + 3 * Mz * nz) * ny) - My);
 
    % 计算位置向量
    r_xyz = r * n_xyz;
    M_xyz = (4*pi*r^4)*M_xyz/(3*mu0);
    % % 使用当前解计算磁场
    % 
    r_mag = norm(r_xyz);
    B_calculated = (mu0/(4*pi)) * (3*(dot(M_xyz,r_xyz)*r_xyz)/(r_mag^5) - M_xyz/(r_mag^3));

    

    
    % 计算与原始磁场的误差
    B_original = B;
    error = norm(B_calculated - B_original);
    
    % 存储结果
    errors(combo) = error;
    all_M_xyz(:,combo) = M_xyz;
    all_r_xyz(:,combo) = r_xyz;
end

% 找到误差最小的组合
[min_error, best_combo] = min(errors);

% 使用最佳组合的结果
M_xyz_final = all_M_xyz(:,best_combo);
r_xyz_final = all_r_xyz(:,best_combo);
%偶极子坐标转换
r_diople = r_obs - r_xyz_final;
% 输出最终结果
disp('最佳解的误差：');
disp(min_error);
disp('最终磁矩向量 M_xyz：');
disp(M_xyz_final);
disp('最终位置向量 r_xyz：');
disp(r_xyz_final);

% 验证最终结果
r_mag_final = norm(r_xyz_final);
B_final = (mu0/(4*pi)) * (3*(dot(M_xyz_final,r_xyz_final)*r_xyz_final)/(r_mag_final^5) - M_xyz_final/(r_mag_final^3));
disp('计算得到的磁场：');
disp(B_final);
disp('原始磁场：');
disp([Bx; By; Bz]);
disp('相对误差（%）：');
disp(100 * norm(B_final - [Bx; By; Bz]) / norm([Bx; By; Bz]));