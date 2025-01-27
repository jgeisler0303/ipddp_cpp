function f = step(in1)
%STEP
%    F = STEP(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    24-Jun-2024 12:42:45

u1 = in1(5,:);
u2 = in1(6,:);
x_1 = in1(1,:);
x_2 = in1(2,:);
x_3 = in1(3,:);
x_4 = in1(4,:);
t2 = cos(u1);
t3 = sin(u1);
t4 = x_4.^2;
t5 = t3.^2;
t6 = t2.*x_4.*(3.0./1.0e+2);
t7 = t4.*t5.*9.0e-4;
t8 = -t7;
t9 = t8+4.0;
t10 = sqrt(t9);
t11 = -t10;
t12 = t6+t11+2.0;
f = [x_1+t12.*cos(x_3);x_2+t12.*sin(x_3);x_3+asin(t3.*x_4.*(3.0./2.0e+2));u2.*(3.0./1.0e+2)+x_4];
