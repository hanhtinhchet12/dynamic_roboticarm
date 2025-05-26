% Các tham số
clear all;clc;
syms theta1 theta2 theta3 real
syms theta1_d theta2_d theta3_d real
syms J1 J2 J3 l_1 l_2 l_3 m1 m2 m3
theta = [theta1; theta2; theta3];
theta_d = [theta1_d; theta2_d; theta3_d];


% Định nghĩa M(θ)
M11 = (J1 + J2 + J3 + (l_1^2*m1)/4 + l_1^2*m2 + l_1^2*m3 + (l_2^2*m2)/4 + l_2^2*m3 + (l_3^2*m3)/4 + l_1*l_3*m3*cos(theta2 + theta3) + l_1*l_2*m2*cos(theta2) + 2*l_1*l_2*m3*cos(theta2) + l_2*l_3*m3*cos(theta3));

M12 = (J2 + J3 + (l_2^2*m2)/4 + l_2^2*m3 + (l_3^2*m3)/4 + (l_1*l_3*m3*cos(theta2 + theta3))/2 + (l_1*l_2*m2*cos(theta2))/2 + l_1*l_2*m3*cos(theta2) + l_2*l_3*m3*cos(theta3));

M13 = (J3 + (l_3^2*m3)/4 + (l_1*l_3*m3*cos(theta2 + theta3))/2 + (l_2*l_3*m3*cos(theta3))/2);

M21 = M12;
M22 = (J2 + J3 + (l_2^2*m2)/4 + l_2^2*m3 + (l_3^2*m3)/4 + l_2*l_3*m3*cos(theta3));
M23 = ((m3*l_3^2)/4 + (l_2*m3*cos(theta3)*l_3)/2 + J3);

M31 = M13;
M32 = M23;
M33 = ((m3*l_3^2)/4 + J3);

M = [M11 M12 M13;
     M21 M22 M23;
     M31 M32 M33];

n = 3;
C = sym(zeros(n, n));
for i = 1:n
    for j = 1:n
        cij = 0;
        for k = 1:n
            cij = cij + 0.5 * ( diff(M(i,j), theta(k)) + ...
                                diff(M(i,k), theta(j)) - ...
                                diff(M(j,k), theta(i)) ) * theta_d(k);
        end
        C(i,j) = simplify(cij);
    end
end

% Tính vector C(theta, theta_dot) * theta_dot
Cqdot = simplify(C * theta_d);

disp('Ma trận Coriolis C(q, q̇):');
disp(C);

disp('Vector Coriolis và ly tâm C(q, q̇)*q̇:');
disp(Cqdot);
