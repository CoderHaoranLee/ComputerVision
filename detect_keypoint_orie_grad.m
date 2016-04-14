function [ keypoint_location,keypoint_orie,keypoint_grad,pyramid_L,pyramid_D,extrema_pos ] = detect_keypoint_orie_grad( img )
%               img: ����Ҷ�ͼ��
% keypoint_location: ����ؼ���λ�ã���matlab����ϵ��λ�ã���ֱ������ͼ���ϻ���
%     keypoint_orie: ����ؼ����������
%     keypoint_grad: ����ؼ�����ݶ�ֵ
%         pyramid_L: ���������LoG������
%         pyramid_D: ���������DoG������
%       extrema_pos: �����ֵ���ڽ������е�λ��[octave,s,m,n]
% ˵���� �ð汾��Ϊ��һ�汾��ͬ־�ǿ���ʱ�øð汾������������������������δ�����
%      1.��ֵ�õ���octave�м���ļ�ֵ��λ��ƫ��Ƚϴ�����δ������뼫ֵ��
%      2.ȥ����������ĵļ�ֵ���ʣ�¼�ֵ��λ��ƫ��Ƚ�С����δ����λ�������Լ������ؼ���λ
%      3.�ڼ���������ʱ��δ��������ֵ��80%����ֵ��Ϊ��������
%      4.������λ���е�����ò��
img = im2double(img);

% origional image gaussian blur parameter
sigma_orig = 0.5;
gauss_size = floor(sigma_orig*6+1); % 3 sigma principle down intergrate
gauss_func = fspecial('gaussian',[gauss_size,gauss_size],sigma_orig);

% gaussian blur and image linear interpolate
img_gauss = conv2(img,gauss_func,'same');
img_linear = imresize(img_gauss,2*size(img));

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
            sigma = sigma_init*sqrt(2)^(o-1 + (s-1)/S);
            gauss_size = floor(sigma*6+1);
            gauss_func = fspecial('gaussian',[gauss_size,gauss_size],sigma);
            octave(:,:,s) = conv2(img_linear,gauss_func,'same');
        end
        pyramid_L{o} = octave;
    else
        for s = 1:S+3
            sigma = sigma_init*sqrt(2)^(o-1 + (s-1)/S);
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
extrema_flag = ones(size(X_org,1),1);
for i = 1:size(extrema_pos,1);
    X = [];
    D = [];
    o = extrema_pos(i,1);
    s = extrema_pos(i,2);
    x = extrema_pos(i,3);
    y = extrema_pos(i,4);
    % ���ز�ֵõ�
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
    % ȥ������ƫ�ƱȽϴ�ĵ㣬��Ҫ�����ڵ�һoctave����ֵ�õ������飬ԭ��δ֪��
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
key_point_descrip = zeros(size(extrema_pos,1),2);
keypoint_orie = [];
keypoint_grad = [];
bin_with = 360/35;
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
                theta_sin = (L_temp(x_i,y_i+1)-L_temp(x_i,y_i-1))/m_i;
                theta_cos = (L_temp(x_i+1,y_i)-L_temp(x_i-1,y_i))/m_i;
                theta = (180/pi)*atan((L_temp(x_i,y_i+1)-L_temp(x_i,y_i-1))/(L_temp(x_i+1,y_i)-L_temp(x_i-1,y_i)));
                % �ѽǶȴ�-180-180ӳ�䵽0-360
                if (theta_sin>0&&theta_cos<0) % �ڶ�����
                    theta = theta + 180;
                elseif (theta_sin<0&&theta_cos<0) % ��������
                    theta = theta + 180;
                elseif (theta_sin<0&&theta_cos>0) % ��������
                    theta = theta + 360;
                end
                if(theta_sin==1)
                    theta = 90;
                elseif(theta_sin == -1)
                    theta = 270;
                end
                bin_i = floor(theta/bin_with);
                delta_theta = theta - bin_i*bin_with;
                bin_i = bin_i + 1;
                delta_add = L_temp(x_i,y_i)*(1/(sqrt(2*pi)*sigma))*exp(((x_i-x)^2+(y_i-y)^2)/(2*sigma^2));
                if(delta_theta<bin_with/2)
                    key_point_hist(i,bin_i) = key_point_hist(i,bin_i) + delta_add;
                else
                    bin_i = bin_i + 1;
                    key_point_hist(i,bin_i) = key_point_hist(i,bin_i) + delta_add;
                end
            end
        end
        [~,id] = max(key_point_hist(i,:));
        theta = bin_with*(id-1);
        keypoint_orie = [keypoint_orie;theta];
        m = sqrt((L_temp(x+1,y_i) - L_temp(x-1,y))^2+(L_temp(x,y+1) - L_temp(x,y-1))^2);
        keypoint_grad = [keypoint_grad;m];
        key_point_descrip(i,1) = m*sin(theta);
        key_point_descrip(i,2) = m*cos(theta);
    end
end
key_point_descrip = key_point_descrip(edge_flag==1,:);
key_point_hist = key_point_hist(edge_flag==1,:);
% =========================================================================
% ����ȥ��ͼ����ش��ĵ�֮����������Լ��õ㴦���ݶ�ʸ��
extrema_pos_org = extrema_pos_org(edge_flag==1,:);
figure(5)
imshow(img)
hold on
plot(extrema_pos_org(:,2),extrema_pos_org(:,1),'r+')
hold on
quiver(extrema_pos_org(:,2),extrema_pos_org(:,1),key_point_descrip(:,1),key_point_descrip(:,2),2,'color', 'y')
% =========================================================================
keypoint_location = [extrema_pos_org(:,2),extrema_pos_org(:,1)];
extrema_pos = extrema_pos(edge_flag==1,:);
end

