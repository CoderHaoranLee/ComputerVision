% match keypoints
clear
close all
clc
load keypoint_descrip_1
load keypoint_descrip_2
load keypoint_location1
load keypoint_location2
[m1,~] = size(keypoint_descrip_1);
[m2,~] = size(keypoint_descrip_2);
if m1>m2
    template = keypoint_descrip_1;
    goal = keypoint_descrip_2;
    m_temp = m1;
    m_g = m2;
    loca_temp = keypoint_location1;
    loca_g = keypoint_location2;
else
    template = keypoint_descrip_2;
    goal = keypoint_descrip_1;
    m_temp = m2;
    m_g = m1;
    loca_temp = keypoint_location2;
    loca_g = keypoint_location1;
end
match_result = zeros(m_g,1);
radio_th = 0.6;
kd_tree = kdtree_build( template );
for i = 1:m_g
    keypoint = goal(i,:);
%     keypoint_M = repmat(keypoint,m_temp,1);
%     delta_M = keypoint_M - template;
%     distance = sqrt(diag(delta_M*delta_M'));
    idxs = kdtree_k_nearest_neighbors(kd_tree,keypoint,2);
    distance_1 = sqrt((keypoint - template(idxs(1),:))* (keypoint - template(idxs(1),:))');
    distance_2 = sqrt((keypoint - template(idxs(2),:))* (keypoint - template(idxs(2),:))');
    [distance_sort,id] = sort([distance_1,distance_2]);
    radio = distance_sort(1)/distance_sort(2);
    if radio<radio_th
        match_result(i) = idxs(id(1));
    end
end
img1 = imread('1.JPG');
img1 = rgb2gray(img1);
img1 = im2double(img1);
img2 = imread('2.JPG');
img2 = rgb2gray(img2);
img2 = im2double(img2);
img = cat(2,img1,img2);
[m,n] = size(img1);
figure
imshow(img)
hold on
loca_g(match_result==0,:) = [];
match_result(match_result==0) = [];
loca_temp = loca_temp(match_result,:);
loca_temp(:,1) = loca_temp(:,1) + repmat(n,size(loca_temp,1),1);
plot(loca_g(:,1),loca_g(:,2),'y+','linewidth',2)
plot(loca_temp(:,1),loca_temp(:,2),'y+','linewidth',2)
for i = 1:size(loca_g,1)
    plot([loca_g(i,1),loca_temp(i,1)],[loca_g(i,2),loca_temp(i,2)],'c-')
end
disp('The result number is')
length(match_result)
