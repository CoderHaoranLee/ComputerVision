function [ output_args ] = harris_detect( img_name )

img = imread(img_name);
img = rgb2gray(img);
img = im2double(img);
sobel_x = [-1 0 1;-2 0 2;-1 0 1];
sobel_y = [1 2 1;0 0 0;-1 -2 -1];
img_grad_x = conv2(img,sobel_x,'same');
img_grad_y = conv2(img,sobel_y,'same');
h = fspecial('gaussian',7,2); % guassian filter, size:3X3, sigma = 2
Ix2 = img_grad_x.^2;
Iy2 = img_grad_y.^2;
Ixy = img_grad_x.*img_grad_y;
M11 = conv2(Ix2,h,'same');
M22 = conv2(Iy2,h,'same');
M12 = conv2(Ixy,h,'same');
R = zeros(size(img));
k = 0.06;
for i = 1:size(img,1)
    for j = 1:size(img,2)
        M = [M11(i,j),M12(i,j);M12(i,j),M11(i,j)];
        R(i,j) = det(M) - k*(trace(M));
    end
end
% find out the maximume in 8 neighborholds
R_max = zeros(size(R));
mean_R = 0.05*max(max(R));
for i = 2:size(img,1)-1
    for j = 2:size(img,2)-1
        if R(i,j)>mean_R
            if (R(i,j)>=R(i-1,j-1))&&(R(i,j)>=R(i,j-1))&&(R(i,j)>=R(i+1,j-1))...
                    &&(R(i,j)>=R(i-1,j))&&(R(i,j)>=R(i+1,j))&&(R(i,j)>=R(i-1,j+1))...
                    &&(R(i,j)>=R(i,j+1))&&(R(i,j)>=R(i+1,j+1))
                R_max(i,j) = 1;
            end
        end
    end
end
R_m = R(R_max==1);
[m,n] = find(R_max==1);
figure
imshow(img)
hold on
plot(n,m,'ro')
% imagesc(R)
% pause(10);


end

