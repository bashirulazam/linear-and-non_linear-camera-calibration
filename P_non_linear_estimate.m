function P = P_non_linear_estimate(P_0,points_3D,points_2D)


P_new = P_0;
p1_new = P_new(1,1:3)';
p2_new = P_new(2,1:3)';
p3_new = P_new(3,1:3)';

p14_new = P_new(1,4);
p24_new = P_new(2,4);
p34_new = P_new(3,4);

N = size(points_3D,2);
i = 1; 

ita = 1e-9;
ita1 = 1e-3;
lambda1 = 1e-04;
lambda2 = 1e-04;

count = 1;
while(count < 1000)
  P = P_new;  
  p1 = p1_new;
  p2 = p2_new;
  p3 = p3_new;
  p14 = p14_new;
  p24 = p24_new;
  p34 = p34_new;
    
    
    for i = 1:1:N
        M(:,i) = points_3D(:,i);
        ch(i) = (M(:,i)'*p1 + p14)./(M(:,i)'*p3 + p34);
        rh(i) = (M(:,i)'*p2 + p24)./(M(:,i)'*p3 + p34);
        c(i) = points_2D(1,i);
        r(i) = points_2D(2,i);
        V(2*i - 1,:) = c(i);
        V(2*i,:) = r(i);
        Vh(2*i - 1,:) = ch(i);
        Vh(2*i,:) = rh(i);
        %for p3
        delc_delp3(:,i) = -((M(:,i)'*p1 + p14)*M(:,i)'*eye(3))./(M(:,i)'*p3 + p34)^2;
        delr_delp3(:,i) = -((M(:,i)'*p2 + p24)*M(:,i)'*eye(3))./(M(:,i)'*p3 + p34)^2;
        delvh_delp3(:,2*i - 1) = delc_delp3(:,i);
        delvh_delp3(:,2*i )    = delr_delp3(:,i);

        %for p2
        delc_delp2(:,i) = zeros(1,3);
        delr_delp2(:,i) = (M(:,i)'*eye(3))./(M(:,i)'*p3 + p34);
        delvh_delp2(:,2*i - 1) = delc_delp2(:,i);
        delvh_delp2(:,2*i) = delr_delp2(:,i);

        %for p1
        delc_delp1(:,i) = (M(:,i)'*eye(3))./(M(:,i)'*p3 + p34);
        delr_delp1(:,i) = zeros(1,3);
        delvh_delp1(:,2*i - 1) = delc_delp1(:,i);
        delvh_delp1(:,2*i) = delr_delp1(:,i);

        %for p14
        delc_delp14(:,i) = 1./(M(:,i)'*p3 + p34);
        delr_delp14(:,i) = 0 ;
        delvh_delp14(:,2*i - 1) = delc_delp14(:,i);
        delvh_delp14(:,2*i) =  delr_delp14(:,i);

        %for p24
        delc_delp24(:,i) = 0 ;
        delr_delp24(:,i) = 1./(M(:,i)'*p3 + p34);
        delvh_delp24(:,2*i - 1) = delc_delp24(:,i);
        delvh_delp24(:,2*i) =  delr_delp24(:,i);

        %for p34 
        delc_delp34(:,i) = -(M(:,i)'*p1 + p14)./(M(:,i)'*p3 + p34)^2;
        delr_delp34(:,i) = -(M(:,i)'*p2 + p24)./(M(:,i)'*p3 + p34)^2;
        delvh_delp34(:,2*i - 1) = delc_delp34(:,i);
        delvh_delp34(:,2*i) = delr_delp34(:,i);     
    end

    %calculation of  nabla_p3
    I = eye(3);
    p1_cross_delp3_delp3 = [cross(p1,I(:,1)) cross(p1,I(:,2)) cross(p1,I(:,3))];
    lambda2_p3_1st = [p1_cross_delp3_delp3(:,1)'*cross(p2,p3); p1_cross_delp3_delp3(:,2)'*cross(p2,p3); p1_cross_delp3_delp3(:,3)'*cross(p2,p3)] ; 
    p2_cross_delp3_delp3 = [cross(p2,I(:,1)) cross(p2,I(:,2)) cross(p2,I(:,3))];
    lambda2_p3_2nd = [p2_cross_delp3_delp3(:,1)'*cross(p1,p3); p2_cross_delp3_delp3(:,2)'*cross(p1,p3); p2_cross_delp3_delp3(:,3)'*cross(p1,p3)] ;
    lambda2_p3 = lambda2_p3_1st + lambda2_p3_2nd;
    nabla_p3 = 2*delvh_delp3*(Vh-V) + 2*lambda1*p3 + lambda2*lambda2_p3;

    %calculation of nabla_p2
    p3_cross_delp2_delp2 = [cross(p3,I(:,1)) cross(p3,I(:,2)) cross(p3,I(:,3))];
    lambda2_p2 = [p3_cross_delp2_delp2(:,1)'*cross(p1,p3) ; p3_cross_delp2_delp2(:,2)'*cross(p1,p3) ;p3_cross_delp2_delp2(:,3)'*cross(p1,p3)];
    nabla_p2 = 2*delvh_delp2*(Vh-V) + lambda2*lambda2_p2;

    %calculation of nabla_p1
    p3_cross_delp1_delp1 = [cross(p3,I(:,1)) cross(p3,I(:,2)) cross(p3,I(:,3))];
    lambda2_p1 = [p3_cross_delp1_delp1(:,1)'*cross(p2,p3) ; p3_cross_delp1_delp1(:,2)'*cross(p2,p3) ; p3_cross_delp1_delp1(:,2)'*cross(p2,p3)];
    nabla_p1 = 2*delvh_delp1*(Vh-V) + lambda2*lambda2_p1;

    %calculation of nabla_p14, nabla_p23, nabla_p34 
    nabla_p14 = 2*delvh_delp14*(Vh - V);
    nabla_p24 = 2*delvh_delp24*(Vh - V);
    nabla_p34 = 2*delvh_delp34*(Vh - V);

    p3_new = p3 + ita*-nabla_p3;
    p2_new = p2 + ita1*-nabla_p2;
    p1_new = p1 + ita1*-nabla_p1;
    p14_new = p14 + ita1*-nabla_p14;
    p24_new = p24 + ita1*-nabla_p24;
    p34_new = p34 + ita1*-nabla_p34;
    
    P_new(1,1:3) = p1_new';
    P_new(2,1:3) = p2_new';
    P_new(3,1:3) = p3_new';
    P_new(:,4) = [ p14_new; p24_new ; p34_new];    
        
    Err(count) = norm(P - P_new,2)/norm(P_new,2);
   
%     if Err(count) < 1e-5 
%         break
%     end
   count = count + 1; 
end
%plot(Err)

