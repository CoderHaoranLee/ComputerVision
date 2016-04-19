% SIFT
clear
close all
clc
img = imread('rice.jpg');
% img = rgb2gray(img);
img = im2double(img);

% origional image gaussian blur parameter
sigma_orig = 0.5;
gauss_size = floor(sigma_orig*6+1); % 3 sigma principle down intergrate
gauss_func = fspecial('gaussian',[gauss_size,gauss_size],sigma_orig);

% gaussian blur and image linear interpolate
img_gauss = conv2(img,gauss_func,'same');
% img_linear = imresize(img_gauss,2*size(img));
img_linear = doubleSize(img_gauss);
%% construct scale space

% scale space parameter and construct LoG pyramid
S = 3;
sigma_init = 1.6;
octave_num =4;  % the size of image on the top of pyramid is 12*16
% pyramid_L = {};
[m,n] = size(img_linear);
for o = 1:octave_num
    octave = zeros(m/2^(o-1),n/2^(o-1),S+3);
    if (o==1)
        for s = 1:S+3
            sigma = sigma_init*2^(o-1 + (s-1)/S);
            gauss_size = floor(sigma*6+1);
            gauss_func = fspecial('gaussian',[gauss_size,gauss_size],sigma);
            octave(:,:,s) = conv2(img_linear,gauss_func,'same');
        end
        pyramid_L{o} = octave;
    else
        for s = 1:S+3
            sigma = sigma_init*2^(o-1 + (s-1)/S);
            gauss_size = floor(sigma*6+1);
            gauss_func = fspecial('gaussian',[gauss_size,gauss_size],sigma);
            img_temp = pyramid_L{o-1}(:,:,S+1);
            img_temp = imresize(img_temp,size(img_temp)/2);
            octave(:,:,s) = conv2(img_temp,gauss_func,'same');
        end
        pyramid_L{o} = octave;
    end
end

% construct DoG pyramid
for o = 1:octave_num
    [m,n,r] = size(pyramid_L{o});
    octave = zeros(m,n,r-1);
    for s = 1:S+2
        img_up = pyramid_L{o}(:,:,s+1);
        img_down = pyramid_L{o}(:,:,s);
        octave(:,:,s) = img_up-img_down;
    end
    pyramid_D{o} = octave;
end
%% plot
plotDoG(pyramid_D)

%% detect extrema
extrema_pos = [];
for o = 1:octave_num
    for s = 2:S+1
        temp_down = pyramid_D{o}(:,:,s-1);
        temp = pyramid_D{o}(:,:,s);
        temp_up = pyramid_D{o}(:,:,s+1);
        [m,n] = size(temp);
        for x = 2:m-1
            for y = 2:n-1
                down_max = max(max(temp_down(x-1:x+1,y-1:y+1)));
                down_min = min(min(temp_down(x-1:x+1,y-1:y+1)));
                up_max = max(max(temp_up(x-1:x+1,y-1:y+1)));
                up_min = min(min(temp_up(x-1:x+1,y-1:y+1)));
                temp_left_max = max(temp(x-1,y-1:y+1));
                temp_left_min = min(temp(x-1,y-1:y+1));
                temp_right_max = max(temp(x+1,y-1:y+1));
                temp_right_min = min(temp(x+1,y-1:y+1));
                if(temp(x,y)>temp(x,y+1)&&temp(x,y)>temp(x,y-1)&&temp(x,y)>temp_left_max...
                        &&temp(x,y)>temp_right_max&&temp(x,y)>down_max&&temp(x,y)&&up_max)
                    extrema_pos = [extrema_pos;o,s,x,y];
                end
                if(temp(x,y)<temp(x,y+1)&&temp(x,y)<temp(x,y-1)&&temp(x,y)<temp_left_min...
                        &&temp(x,y)<temp_right_min&&temp(x,y)<down_min&&temp(x,y)&&up_min)
                    extrema_pos = [extrema_pos;o,s,x,y];
                end
            end
        end
    end
end
% =========================================================================
% ����������ļ�ֵ��
extrema_pos_org = [extrema_pos(:,end-1).*(2.^(extrema_pos(:,1)-2)),extrema_pos(:,end).*(2.^(extrema_pos(:,1)-2))];
figure(2)
imshow(img)
hold on
plot(extrema_pos_org(:,2),extrema_pos_org(:,1),'r+')
% =========================================================================

%%  fit the quadratic function
Sigma = sigma_init*2.^(extrema_pos(:,1)-1 + (extrema_pos(:,2)-1)/S);
X_org = [extrema_pos(:,end-1),extrema_pos(:,end),Sigma];
delta_X = [];
delta_X_pixel = [];

% ���ֹؼ��㾫ȷ��λ�ķ���

extrema_flag = ones(size(X_org,1),1);
for i = 1:size(extrema_pos,1);
    X = [];
    D = [];
    o = extrema_pos(i,1);
    s = extrema_pos(i,2);
    x = extrema_pos(i,3);
    y = extrema_pos(i,4);
    
    % 1.����ÿ�������㣬Ҫ�ҳ���26�����ֵ�����Ԫ���κ�����ȥ���ԱȶȱȽϵ͵ĵ�
    %======================================================================
    %     % 3 layer
    %     for s_i = s-1:s+1
    %         for x_i = x-1:x+1
    %             for y_i = y-1:y+1
    %                 sigma = sigma_init*2^(o-1 + (s_i-1)/S);
    %                 X_temp = [x_i,y_i,sigma]- X_org(i,:);
    %                 X = [X;X_temp];
    %                 D = [D;pyramid_D{o}(x_i,y_i,s_i)];
    %             end
    %         end
    %     end
    %     x1 = X(:,1);
    %     x2 = X(:,2);
    %     x3 = X(:,3);
    %     x4 = x1.^2;
    %     x5 = x1.*x2;
    %     x6 = x1.*x3;
    %     x7 = x2.^2;
    %     x8 = x2.*x3;
    %     x9 = x3.^2;
    %     X_fit = [ones(size(x1)),x1,x2,x3,x4,x5,x6,x7,x8,x9];
    %     b = pinv(X_fit'*X_fit)*X_fit'*D;
    %     d = b(1);
    %     d1 = b(2);
    %     d2 = b(3);
    %     d3 = b(4);
    %     d11 = b(5);
    %     d12 = b(6);
    %     d13 = b(7);
    %     d22 = b(8);
    %     d23 = b(9);
    %     d33 = b(10);
    %     D_g = [d1;d2;d3];
    %     D_h = 2*[d11,d12,d13;d12,d22,d23;d13,d23,d33];
    %     delta_x = -pinv(D_h)*D_g;
    % %     if abs(d+0.5*D_g'*delta_x)<0.04
    % %         extrema_flag(i) = 0;
    % %     end
    %     delta_X = [delta_X;delta_x'];
    %======================================================================
    % 2.���ز�ֵõ�
    % ���޲���󵼲���
    deriv_scale = 1/2;
    second_deriv_scale = 1;
    cross_deriv_scale = 1/4;
    D_temp = pyramid_D{o};
    dx = (D_temp(x+1,y,s) - D_temp(x-1,y,s))*deriv_scale;
    dy = (D_temp(x,y+1,s) - D_temp(x,y-1,s))*deriv_scale;
    ds = (D_temp(x,y,s+1) - D_temp(x,y,s-1))*deriv_scale;
    d_temp = D_temp(x,y,s)*2;
    dxx = (D_temp(x+1,y,s) + D_temp(x-1,y,s) - d_temp)*second_deriv_scale;
    dyy = (D_temp(x,y+1,s) + D_temp(x,y-1,s) - d_temp)*second_deriv_scale;
    dss = (D_temp(x,y,s+1) + D_temp(x,y,s-1) - d_temp)*second_deriv_scale;
    dxy = (D_temp(x+1,y+1,s) + D_temp(x-1,y-1,s) - D_temp(x+1,y-1,s) - D_temp(x-1,y+1,s))*cross_deriv_scale;
    dxs = (D_temp(x+1,y,s+1) + D_temp(x-1,y,s-1) - D_temp(x+1,y,s-1) - D_temp(x-1,y,s+1))*cross_deriv_scale;
    dys = (D_temp(x,y+1,s+1) + D_temp(x,y-1,s-1) - D_temp(x,y+1,s-1) - D_temp(x,y-1,s+1))*cross_deriv_scale;
    D_g_pixel = [dx;dy;ds];
    D_h_pixel = [dxx,dxy,dxs;dxy,dyy,dys;dxs,dys,dss];
    delta_x_pixel = -pinv(D_h_pixel)*D_g_pixel;
    delta_X_pixel = [delta_X_pixel;delta_x_pixel'];
    if abs(d_temp+0.5*D_g_pixel'*delta_x_pixel)<0.04
        extrema_flag(i) = 0;
    end
%     % ȥ������ƫ�ƱȽϴ�ĵ㣬��Ҫ�����ڵ�һoctave����ֵ�õ������飬ԭ��δ֪��
    if sum(abs(delta_x_pixel))>5
        extrema_flag(i) = 0;
    end
end

X_org = [X_org,delta_X,delta_X_pixel];
X_detect = X_org(extrema_flag==1,:);
extrema_pos = extrema_pos(extrema_flag==1,:);
% =========================================================================
% ����ȥ���ԱȶȱȽϵ��Լ�ƫ�ƱȽϴ�ĵ�
extrema_pos_org = [extrema_pos(:,end-1).*(2.^(extrema_pos(:,1)-2)),extrema_pos(:,end).*(2.^(extrema_pos(:,1)-2))];
figure(3)
imshow(img)
hold on
plot(extrema_pos_org(:,2),extrema_pos_org(:,1),'r+')
% =========================================================================
%% ȥ����ԵЧӦ
r = 10;
extrema_flag = ones(size(extrema_pos,1),1);
for i = 1:size(extrema_pos,1)
    o = extrema_pos(i,1);
    s = extrema_pos(i,2);
    x = extrema_pos(i,3);
    y = extrema_pos(i,4);
    % ��Hessian����
    second_deriv_scale = 1;
    cross_deriv_scale = 1/4;
    D_temp = pyramid_D{o};
    d_temp = D_temp(x,y,s)*2;
    dxx = (D_temp(x+1,y,s) + D_temp(x-1,y,s) - d_temp)*second_deriv_scale;
    dyy = (D_temp(x,y+1,s) + D_temp(x,y-1,s) - d_temp)*second_deriv_scale;
    dxy = (D_temp(x+1,y+1,s) + D_temp(x-1,y-1,s) - D_temp(x+1,y-1,s) - D_temp(x-1,y+1,s))*cross_deriv_scale;
    H_trace = dxx + dyy;
    H_det = dxx*dyy - dxy^2;
    if (H_trace^2/H_det)>((r+1)^2/r)
        extrema_flag(i) = 0;
    end
end
extrema_pos = extrema_pos(extrema_flag==1,:);
% =========================================================================
% ����ȥ����Ե�ĵ�
extrema_pos_org = [extrema_pos(:,end-1).*(2.^(extrema_pos(:,1)-2)),extrema_pos(:,end).*(2.^(extrema_pos(:,1)-2))];
figure(4)
imshow(img)
hold on
plot(extrema_pos_org(:,2),extrema_pos_org(:,1),'r+')
% =========================================================================
%% ����ؼ���������
edge_flag = ones(size(extrema_pos,1),1);
key_point_hist = zeros(size(extrema_pos,1),36);
key_point_grad = zeros(size(extrema_pos,1),2);
key_point_orie = zeros(size(extrema_pos,1),1);
bin_width = 360/35;
for i = 1:size(extrema_pos,1)
    o = extrema_pos(i,1);
    s = extrema_pos(i,2);
    x = extrema_pos(i,3);
    y = extrema_pos(i,4);
    sigma = sigma_init*2^(o-1 + (s-1)/S);
    radius = floor(3*1.5*sigma);
    [m,n] = size(pyramid_L{o}(:,:,s));
    L_temp = pyramid_L{o}(:,:,s);
    if((x-radius)<=1||x+radius>m-1||(y-radius)<=1||y+radius>n-1) % ȥ��ͼ����ش��ĵ�
        edge_flag(i) = 0;
        continue;
    else
        for x_i = (x-radius):(x+radius)
            for y_i = (y-radius):(y+radius)
                % ����õ��ݶȺͷ���
                m_i = sqrt((L_temp(x_i+1,y_i) - L_temp(x_i-1,y_i))^2+(L_temp(x_i,y_i+1) - L_temp(x_i,y_i-1))^2);
                theta_cos = (L_temp(x_i,y_i+1)-L_temp(x_i,y_i-1))/m_i;
                theta_sin = (L_temp(x_i+1,y_i)-L_temp(x_i-1,y_i))/m_i;
                theta = (180/pi)*atan(theta_sin/theta_cos);
                % �ѽǶȴ�-180-180ӳ�䵽0-360
                if (theta_cos>0&&theta_sin<0) % �ڶ�����
                    theta = theta + 180;
                elseif (theta_cos<0&&theta_sin<0) % ��������
                    theta = theta + 180;
                elseif (theta_cos<0&&theta_sin>0) % ��������
                    theta = theta + 360;
                end
                if(theta_cos==1)
                    theta = 90;
                elseif(theta_cos == -1)
                    theta = 270;
                end
                bin_i = floor(theta/bin_width);
                delta_theta = theta - bin_i*bin_width;
                bin_i = bin_i + 1;
                delta_add = m_i*(1/(sqrt(2*pi)*sigma))*exp(-((x_i-x)^2+(y_i-y)^2)/(2*sigma^2));
                if(delta_theta<bin_width/2)
                    key_point_hist(i,bin_i) = key_point_hist(i,bin_i) + delta_add;
                else
                    bin_i = bin_i + 1;
                    key_point_hist(i,bin_i) = key_point_hist(i,bin_i) + delta_add;
                end
            end
        end
        [max_v,id] = max(key_point_hist(i,:));
        theta = bin_width*(id-1);
        key_point_orie(i) = theta;
        m = sqrt((L_temp(x+1,y) - L_temp(x-1,y))^2+(L_temp(x,y+1) - L_temp(x,y-1))^2);
        key_point_grad(i,1) = m*sin(theta);
        key_point_grad(i,2) = m*cos(theta);
    end
end
key_point_grad = key_point_grad(edge_flag==1,:);
key_point_hist = key_point_hist(edge_flag==1,:);
key_point_orie = key_point_orie(edge_flag==1,:);
% =========================================================================
% ����ȥ��ͼ����ش��ĵ�֮����������Լ��õ㴦���ݶ�ʸ��
extrema_pos_org = extrema_pos_org(edge_flag==1,:);
figure(5)
imshow(img)
hold on
plot(extrema_pos_org(:,2),extrema_pos_org(:,1),'r+')
hold on
quiver(extrema_pos_org(:,2),extrema_pos_org(:,1),key_point_grad(:,1),key_point_grad(:,2),2,'color', 'y')
% =========================================================================
%% �����ؼ���������
d = 4;
bins = 8;
bin_width = 360/(bins-1);
key_point_descrip = zeros(size(extrema_pos_org,1),d*d*bins);
for i = 1:size(extrema_pos_org,1)
    o = extrema_pos(i,1);
    s = extrema_pos(i,2);
    x = extrema_pos(i,3);
    y = extrema_pos(i,4);
    sigma = sigma_init*2^(o-1 + (s-1)/S);
    w = floor(3*sigma);
    radius_d = (d-1)*w/2;
    descriptor = zeros(bins,d*d);
    descriptor_normal = zeros(d*d*bins,1);
    L_temp = pyramid_L{o}(:,:,s);
    [m,n] = size(L_temp);
    theta_org = key_point_orie(i)*pi/180;
    RM = [cos(theta_org) -sin(theta_org);sin(theta_org) cos(theta_org)];
    k = 1;
    for x_reg = (x - radius_d):w:(x + radius_d)
        for y_reg = (y - radius_d):w:(y + radius_d)
            x_reg = floor(x_reg);
            y_reg = floor(y_reg);
            hist = zeros(bins,1);
            for x_i = (x_reg-w):(x_reg+w)
                for y_i = (y_reg-w):(y_reg+w)
                    x_o = [x_i - x;y_i - y];
                    x_r = RM*x_o;
                    if((x_i)<=1||x_i>m-1||(y_i)<=1||y_i>n-1) % ȥ��ͼ����ش��ĵ�
                       % hist = hist;
                    else
                        % ����õ��ݶȺͷ���
                        m_i = sqrt((L_temp(x_i+1,y_i) - L_temp(x_i-1,y_i))^2+(L_temp(x_i,y_i+1) - L_temp(x_i,y_i-1))^2);
                        theta_cos = (L_temp(x_i,y_i+1)-L_temp(x_i,y_i-1))/m_i;
                        theta_sin = (L_temp(x_i+1,y_i)-L_temp(x_i-1,y_i))/m_i;
                        theta = (180/pi)*atan((L_temp(x_i,y_i+1)-L_temp(x_i,y_i-1))/(L_temp(x_i+1,y_i)-L_temp(x_i-1,y_i)));
                        % �ѽǶȴ�-180-180ӳ�䵽0-360
                        if (theta_cos>0&&theta_sin<0) % �ڶ�����
                            theta = theta + 180;
                        elseif (theta_cos<0&&theta_sin<0) % ��������
                            theta = theta + 180;
                        elseif (theta_cos<0&&theta_sin>0) % ��������
                            theta = theta + 360;
                        end
                        if(theta_cos==1)
                            theta = 90;
                        elseif(theta_cos == -1)
                            theta = 270;
                        end
                        bin_i = floor(theta/bin_width);
                        delta_theta = theta - bin_i*bin_width;
                        bin_i = bin_i + 1;  % ��������ֱ��ͼλ��
                        do = delta_theta/bin_width;
                        dr = abs(y_i - y_reg)/bin_width;
                        dc = abs(x_i - x_reg)/bin_width;
                        delta_add = do*dr*dc*m_i*(1/(sqrt(2*pi)*0.5*d*w))*exp(-(x_r'*x_r)/(2*(0.5*d*w)^2));
                        delta_add_1 = (1-do)*dr*dc*m_i*(1/(sqrt(2*pi)*0.5*d*w))*exp(-(x_r'*x_r)/(2*(0.5*d*w)^2));
                        hist(bin_i) = hist(bin_i) + delta_add;
                        if(bin_i+1>bins)
                            hist(bin_i-bins) = hist(bin_i-bins) + delta_add;
                        else
                            hist(bin_i+1) = hist(bin_i+1) + delta_add_1;
                        end
                    end
                end
            end
            descriptor(:,k) = hist;
            
            k = k+1;
        end
    end
    if (norm(descriptor(:))~= 0)
        descriptor_normal = descriptor(:)/norm(descriptor(:));
    end
    descriptor_normal(descriptor_normal>0.2) = 0.2;
    if (norm(descriptor(:))~= 0)
        descriptor_normal = descriptor(:)/norm(descriptor(:));
    end
    key_point_descrip(i,:) = descriptor_normal';
end

% ========================================================================
