function [ key_point_descrip ] = construct_descriptors( extrema_pos,key_point_orie,pyramid_L )
%% 构建关键点描述子
S = 3;
sigma_init = 1.6;
d = 4;
bins = 8;
bin_width = 360/(bins-1);
key_point_descrip = zeros(size(extrema_pos,1),d*d*bins);
for i = 1:size(extrema_pos,1)
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
                    if((x_i)<=1||x_i>m-1||(y_i)<=1||y_i>n-1) % 去掉图像边沿处的点
                       % hist = hist;
                    else
                        % 计算该点梯度和方向
                        m_i = sqrt((L_temp(x_i+1,y_i) - L_temp(x_i-1,y_i))^2+(L_temp(x_i,y_i+1) - L_temp(x_i,y_i-1))^2);
                        theta_cos = (L_temp(x_i,y_i+1)-L_temp(x_i,y_i-1))/m_i;
                        theta_sin = (L_temp(x_i+1,y_i)-L_temp(x_i-1,y_i))/m_i;
                        theta = (180/pi)*atan((L_temp(x_i,y_i+1)-L_temp(x_i,y_i-1))/(L_temp(x_i+1,y_i)-L_temp(x_i-1,y_i)));
                        % 把角度从-180-180映射到0-360
                        if (theta_cos>0&&theta_sin<0) % 第二象限
                            theta = theta + 180;
                        elseif (theta_cos<0&&theta_sin<0) % 第三象限
                            theta = theta + 180;
                        elseif (theta_cos<0&&theta_sin>0) % 第四象限
                            theta = theta + 360;
                        end
                        if(theta_cos==1)
                            theta = 90;
                        elseif(theta_cos == -1)
                            theta = 270;
                        end
                        bin_i = floor(theta/bin_width);
                        delta_theta = theta - bin_i*bin_width;
                        bin_i = bin_i + 1;  % 计算所在直方图位置
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
end

