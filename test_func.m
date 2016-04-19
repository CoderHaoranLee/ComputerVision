clear
close all
clc
img1 = imread('flower1.jpg');
img1 = rgb2gray(img1);
img1 = im2double(img1);
img2 = imread('flower2.jpg');
img2 = rgb2gray(img2);
img2 = im2double(img2);
tic
[ keypoint_location1,keypoint_orie1,keypoint_grad1,pyramid_L1,pyramid_D1,extrema_pos1 ] = detect_keypoint_orie_grad( img1 );
toc
tic
[ keypoint_location2,keypoint_orie2,keypoint_grad2,pyramid_L2,pyramid_D2,extrema_pos2 ] = detect_keypoint_orie_grad( img2 );
toc
tic
keypoint_descrip_1 = construct_descriptors(extrema_pos1,keypoint_orie1,pyramid_L1);
toc
tic
keypoint_descrip_2 = construct_descriptors(extrema_pos2,keypoint_orie2,pyramid_L2);
toc
tic
match_keypoint( keypoint_descrip_1,keypoint_descrip_2,keypoint_location1,keypoint_location2 )
toc