clear
close all
clc
img = imread('rice.jpg');
[ keypoint_location,keypoint_orie,keypoint_grad,pyramid_L,pyramid_D,extrema_pos ] = detect_keypoint_orie_grad( img );
