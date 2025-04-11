% 计算观测平面上任意方向磁矩产生的磁场Bz分量
clear all;
close all;

%% 定义物理常数
mu0 = 4*pi*1e-7;    % 真空磁导率
%% 创建观测平面的网格
no = 100;
x1 = linspace(-20, 20, no);
y1 = linspace(-20, 20, no);
d = 10;            % 偶极子到观测平面的距离
[X, Y] = meshgrid(x1, y1);
Z=d * ones(size(X));

%% 计算磁场分布
Bz = xlsread('F:\A-寒假\磁偶极子空间转换反演\elsarticle\elsarticle\CG投稿文献\代码可用性\text_example_Bz_data.xlsx');
%% 计算一阶梯度
dx = x1(2)-x1(1);
dy = y1(2)-y1(1);
[dBxdx, dBxdy, dBxdz, dBydy, dBydz] = text_Fourier(Bz,dx,dy);
%% 计算并绘制梯度
[dB_dx, dB_dy] = gradient(Bz, x1(2)-x1(1), y1(2)-y1(1));
gradient_amplitude_q = sqrt(dB_dx.^2 + dB_dy.^2);
q_max = max(gradient_amplitude_q(:));
q = gradient_amplitude_q/q_max;
% 找出最大梯度位置
[max_grad_val, max_grad_idx] = max(gradient_amplitude_q(:));
[max_grad_row, max_grad_col] = ind2sub(size(gradient_amplitude_q), max_grad_idx);
max_grad_point = [X(max_grad_row, max_grad_col), Y(max_grad_row, max_grad_col)];
%% 反推
nRows = size(Bz,1);
nCols = size(Bz,2);

%分区
% 找出在该区域内的点的索引
block_size = 10; % 50/5 = 10个点一块
num_blocks = 1;

block = 1;
    % 计算当前块的索引范围
    start_idx =1;
    end_idx = no;
   
    % 提取当前块的数据
    x_block = X(:,start_idx:end_idx);
    y_block = Y(:,start_idx:end_idx);
    z_block = Z(:,start_idx:end_idx);
    B_block = Bz(:,start_idx:end_idx);
    q_block = q(:,start_idx:end_idx);

    % 设置阈值
threshold = 0.7;

%筛选
q_select = q_block > threshold;

    % 获取这些点的坐标值
    x_points = x_block(q_select);
    y_points = y_block(q_select);
    z_points = z_block(q_select);
    B_points = B_block(q_select);
%% 代码画图
[x_estimates,y_estimates,z_estimates,Mx_estimates,My_estimates,Mz_estimates] = text_inverse(x_points,y_points,d,x1,y1,Bz,dBxdx,dBxdy,dBxdz,dBydy,dBydz);
%% 密度函数计算团簇
pra_data = [x_estimates, y_estimates];
% 设置参数，esp为领域大小，minpt为最小数量点
epsilon = 0.43;
minPts = 200;
% 聚类参数
idx = dbscan(pra_data,epsilon,minPts);
% 不包括噪声点（-1）的簇数
num_clusters = length(unique(idx(idx~=-1)));

[top_clusters, top_cluster_data, top_cluster_counts] = find_top_n_clusters(pra_data,idx, 1);

%计算选出点对应的M
num_cluster1_data = length(top_cluster_data);
for i = 1:num_cluster1_data
    % 提取当前簇的数据
    cluster = top_cluster_data{i};
    
    % 提取 x 和 y 坐标
    x_coords{i} = cluster(:, 1);
    y_coords{i} = cluster(:, 2);
end 
%% 计算选出点对应的M
%长方体1
n_1 = length(x_coords{1});
% 初始化cluster1_data_Mx和cluster1_data_My为长度为num_cluster1_data的向量
cluster1_data_Mx = zeros(n_1, 1);
cluster1_data_My = zeros(n_1, 1);
for idx_cluster1 = 1:n_1
    % 找到最接近的x和y索引
    cluster1_data_x = x_coords{1};
    cluster1_data_y = y_coords{1};
    target_cluster1_x = cluster1_data_x(idx_cluster1);
    target_cluster1_y = cluster1_data_y(idx_cluster1);
[~, x_cluster1_idx] = min(abs(x_estimates - target_cluster1_x));
[~, y_cluster1_idy] = min(abs(y_estimates - target_cluster1_y));   
        % 将每次得到的值存储到cluster1_data_Mx和cluster1_data_My的对应位置
    cluster1_data_Mx(idx_cluster1) = Mx_estimates(x_cluster1_idx);
    cluster1_data_My(idx_cluster1) = My_estimates(y_cluster1_idy);
end
% 计算选定数据点的magnitude
magnitude1 = sqrt(cluster1_data_Mx.^2 + cluster1_data_My.^2);
%% 自定义颜色2
custom_colormap = [
   16/255,  70/255,  128/255;   % 深蓝色
   49/255,  124/255, 183/255;   % 蓝色
   109/255, 173/255, 209/255;   % 浅蓝色1
   182/255, 215/255, 232/255;   % 浅蓝色2
   233/255, 241/255, 244/255;   % 最浅蓝色
   251/255, 227/255, 213/255;   % 最浅红色
   246/255, 178/255, 147/255;   % 浅红色1
   220/255, 109/255, 87/255;    % 浅红色2
   183/255, 34/255,  48/255;    % 红色
   109/255, 1/255,   31/255     % 深红色
];
n = 64;  % 新颜色映射中的颜色数量
% 在原有6个颜色之间进行插值，生成n个颜色点
custom_colormap_interp = interp1(1:10, custom_colormap, linspace(1,10,n));

orange = [1, 0.5, 0];  % 或 [255,204,0]/255
%% 绘图
figure(1);
set(gcf,'DefaultLineLineWidth',3);
contourf(X, Y, Bz);
colormap(custom_colormap_interp);
hold on;
plot(x_points, y_points, '.', 'Color',orange, 'MarkerSize', 12);

shading interp;
colorbar;

set(gca,'LineWidth',3,'fontsize',30,'fontname','Times New Roman','FontWeight','bold');
set(gcf,'Position',[20 20 1000 800],'color','w');

figure(2)  % 新建一个图窗
set(gcf,'DefaultLineLineWidth',3);
% 绘制散点图
% 绘制散点图,用小圆球表示数据点
scatter(x_coords{1}, y_coords{1}, 60, magnitude1, 'filled','MarkerEdgeColor', 'k');
colormap(custom_colormap_interp);

% 自定义圆球表面属性
h = findobj(gca,'Type','patch'); 
set(h,'SpecularStrength',0.2,'SpecularExponent',5,'SpecularColorReflectance',0.5);

% 添加光源
light('Position',[1 1 1],'Style','local');
box on;
xlim([-20,20]);
ylim([-20,20]);
grid on;
axis equal;
set(gca,'LineWidth',3,'fontsize',30,'fontname','Times New Roman','FontWeight','bold');
set(gcf,'Position',[20 20 1000 800],'color','w');
