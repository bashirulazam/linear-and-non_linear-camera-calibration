clear all
close all
clc

points_3D_bad = read_points('bad_3dpts.txt');
points_2D = read_points('Left_2Dpoints.txt');
points_3D_good = read_points('3Dpointnew.txt');

P = P_linear_estimate_with_RANSAC(points_3D_good,points_2D);
[W,R,T] = calc_W_R_T(P);
N = size(points_3D_good,2);
for i = 1:N 
    test_point_3D = points_3D_good(:,i);
    test_point_2D = points_2D(:,i);
    lambda(i) = R(3,:)*test_point_3D + T(3);
    est_points_2D = ((P*[test_point_3D;1])./lambda(i));
    est_points_2D = abs(est_points_2D(1:2,1));
    Error(i) = norm(est_points_2D - test_point_2D,2)/norm(test_point_2D,2);
end

plot(Error)
xlabel('Index of points')
ylabel('Error')
title('Error Curve for good dataset with RANSAC');





