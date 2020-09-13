 function [P_refined] = P_linear_estimate_with_RANSAC(points_3D,points_2D)


p = 0.9; % probability that only one subset has no outliers
e = 0.18; % percentage of outliers 
N = size(points_3D,2); % number of data points 
K = 20; % number of points in the subset
threshold_error = 0.1; 
s = round(log(1-p)/log(1-(1-e)^K)); % the number of subsets 
S = {}; % initialized cell array for storing all the sets

for m = 1:s 
    random_set = randperm(N,K); % radomly selected K number
    all_set = 1:N;
    random_points_3D = points_3D(:,random_set); % random K 3D numbers 
    random_points_2D = points_2D(:,random_set); % random K 2D numbers 
    P = compute_P(random_points_3D,random_points_2D); % computation of projection matrix 
    [W,R,T] = calc_W_R_T(P); % computation of intrinsic and extrinsic parameters of P 
    j = 1; 
    inlier = [];
    remaining_set = setdiff(all_set,random_set); % remaining indexes 
    % calculating error for remaining points
    for i = 1:N-K
        test_point_3D = points_3D(:,remaining_set(i));
        test_point_2D = points_2D(:,remaining_set(i));
        lambda(i) = R(3,:)*test_point_3D + T(3);
        test_point_2D_est = ((P*[test_point_3D;1])./lambda(i));
        test_point_2D_est = abs(test_point_2D_est(1:2,1));
        Err(i) = norm(test_point_2D_est - test_point_2D,2)/norm(test_point_2D,2);
        % If error is below threshold, increment the number of inliers 
        if Err(i) < threshold_error
            inlier(j) = remaining_set(i);
            j = j+1;
        end
    end
    inlier = [random_set inlier]; % joining with the initial points
    S{1,m} = inlier; % storing all sets
end

S_length = cellfun('length',S); % getting the lengths of all sets in S 
ind = find(S_length == max(S_length)); 
inlier_set = S{1,ind}; % finding the largest set of S 

inlier_points_3D = points_3D(:,inlier_set);
inlier_points_2D = points_2D(:,inlier_set);
P_refined = compute_P(inlier_points_3D, inlier_points_2D);


