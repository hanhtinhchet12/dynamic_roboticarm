clear all;clc;
syms t theta1(t) theta2(t) theta3(t) theta1_d(t) theta2_d(t) theta3_d(t) theta1_2d(t) theta2_2d(t) theta3_2d(t) l_1 l_2 l_3 m1 m2 m3 g T1 T2 T3 J1 J2 J3

% T1 = input('nhập momen T1 đi : ');
% T2 = input('nhập momen T2 đi : ');
% T3 = input('nhập momen T3 đi : ');
% if isempty(T1)  % Kiểm tra nếu người dùng không nhập gì
%     T1 = 0;    % Gán giá trị mặc định
% end
% if isempty(T2)  % Kiểm tra nếu người dùng không nhập gì
%     T2 = 0;    % Gán giá trị mặc định
% end
% if isempty(T3)  % Kiểm tra nếu người dùng không nhập gì
%     T3 = 0;    % Gán giá trị mặc định
% end

% khop 1
x1 = 0 ;
y1 = 0 ;
z1 =l_1/2;
x1_d = 0;
y1_d = 0;
z1_d = 0;
% khop 2
x2 = l_2*cos(theta2(t))*cos(theta1(t));
y2 = l_2*cos(theta2(t))*sin(theta1(t));
z2 = l_1 + l_2*sin(theta2(t));
x2_d = subs(diff(x2,t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t)]);
y2_d = subs(diff(y2,t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t)]);
z2_d = subs(diff(z2,t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t)]);

% khop 3
x3 = (l_2*cos(theta2(t))+l_3*cos(theta2(t)+theta3(t)))*cos(theta1(t));
y3 = (l_2*cos(theta2(t))+l_3*cos(theta2(t)+theta3(t)))*sin(theta1(t));
z3 = l_1 + l_2*sin(theta2(t))+l_3*sin(theta2(t)+theta3(t));
x3_d = subs(diff(x3,t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t)]);
y3_d = subs(diff(y3,t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t)]);
z3_d = subs(diff(z3,t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t)]);

w1 = theta1_d(t);
w2 = theta2_d(t);
w3 = theta2_d(t)+theta3_d(t);
% dong nang
K1 = 1/2*m1*(x1_d^2 + y1_d^2 + z1_d^2)+1/2*m2*(x2_d^2 + y2_d^2 + z2_d^2)+1/2*m3*(x3_d^2 + y3_d^2 + z3_d^2);
K2 = 1/2*(J1*w1^2+J2*w2^2+J3*w3^2);
K =K1 +K2;
% the nang 
V = g*(m1*1/2*l_1+m2*(l_1 + l_2*sin(theta2(t)))+m3*(l_1 + l_2*sin(theta2(t))+l_3*sin(theta2(t)+theta3(t))));
% ham Euler - Lagrange
L = K - V;

% H1 = dL/dO1
H1 = subs(diff(L,theta1(t)),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t),diff(theta1_d(t), t),diff(theta2_d(t), t),diff(theta3_d(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t)]);
% H2 = dL/dO2
H2 = subs(diff(L,theta2(t)),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t),diff(theta1_d(t), t),diff(theta2_d(t), t),diff(theta3_d(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t)]);
% H3 = dL/dO3
H3 = subs(diff(L,theta3(t)),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t),diff(theta1_d(t), t),diff(theta2_d(t), t),diff(theta3_d(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t)]);


% U1 = (d/dt)(dL/dO1_d)
U1 = subs(diff(diff(L,theta1_d(t)),t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t),diff(theta1_d(t), t),diff(theta2_d(t), t),diff(theta3_d(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t)]);
% U2 = (d/dt)(dL/dO2_d)
U2 = subs(diff(diff(L,theta2_d(t)),t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t),diff(theta1_d(t), t),diff(theta2_d(t), t),diff(theta3_d(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t)]);
% U3 = (d/dt)(dL/dO3_d)
U3 = subs(diff(diff(L,theta3_d(t)),t),[diff(theta1(t), t),diff(theta2(t), t),diff(theta3(t), t),diff(theta1_d(t), t),diff(theta2_d(t), t),diff(theta3_d(t), t)],[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t)]);

T1 = simplify(U1 - H1);
T2 = simplify(U2 - H2);
T3 = simplify(U3 - H3);

T1 = collect(T1,[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t),g]);
T2 = collect(T2,[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t),g]);
T3 = collect(T3,[theta1_d(t),theta2_d(t),theta3_d(t),theta1_2d(t),theta2_2d(t),theta3_2d(t),g]);

% truy tìm hệ số
% hesotheta1_2d = coeffs(T1,theta1_2d(t));
% disp("coeffs T1 = ");disp(hesotheta1_2d(end));
% 
% hesotheta2_2d = coeffs(T1,theta2_2d(t));
% disp("coeffs T1 = ");disp(hesotheta2_2d(end));
% 
% hesotheta3_2d = coeffs(T1,theta3_2d(t));
% disp("coeffs T1 = ");disp(hesotheta3_2d(end));
%disp(K);


% disp("T1 = ");disp(T1);
% disp("T2 = ");disp(T2);
disp("T3 = ");disp(T3);









