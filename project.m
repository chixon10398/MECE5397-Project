clear all
close all
clc

%boundary lengths
a = 0;
b = 2*pi;
%discretize the problem
n = 6; %number of steps to b
h = b/n; %distance between each position
x = 0:h:b;
y = x;
%create the functions
fb = y.*(b-y);
gb = ((b-y).^2).*cos(pi*y/b);
F = zeros(n,n-1); %create vector for F values
U_n = zeros(n,n-1); %U for n
U_n1 = zeros(n,n-1); % U for n+1
%create the elements for F
for j = 1:n
    for i = 1:n-1
        F(j,i) = sin(pi*i/b)*cos((pi/2)*(2*j/b+1));
    end
end

%Boundary conditions
u_x0 = fb;
u_xb = gb;
u_y0 = fb(1)+(x/b).*(gb(1)-fb(1));

for j = 1:n
    for i = 1:n-1
        if j == 1 && i == 1
            U_n1(j,i) = (1/4)*(2*U_n(j+1,i)+U_n(j,i+1))+(h^2/4)*F(j,i)-(1/4)*u_x0(j);
        elseif j ==1 && i < n-1
            U_n1(j,i) = (1/4)*(2*U_n(j+1,i)+U_n(j,i+1)+U_n(j,i-1))+(h^2/4)*F(j,i);
        elseif j == 1 && i == n-1
            U_n1(j,i) = (1/4)*(2*U_n(j+1,i)+U_n(j,i-1))+(h^2/4)*F(j,i)-(1/4)*u_xb(j);
        elseif j > 1 && j < n && i == 1
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n1(j-1,i)+U_n(j,i+1))+(h^2/4)*F(j,i)-(1/4)*u_x0(j);
        elseif j > 1 && j < n && i < n-1 && i >1
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n1(j-1,i)+U_n(j,i+1)+U_n1(j,i-1))+(h^2/4)*F(j,i);
        elseif j > 1 && j < n && i == n-1
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n1(j-1,i)+U_n1(j,i-1))+(h^2/4)*F(j,i)-(1/4)*u_xb(j);
        elseif j == n && i == 1
            U_n1(j,i) = (1/4)*(U_n1(j-1,i)+U_n(j,i+1))+(h^2/4)*F(j,i)-(1/4)*(u_x0(j)+u_y0(i+1));
        elseif j == n && i > 1 && i < n-1
            U_n1(j,i) = (1/4)*(U_n1(j-1,i)+U_n(j,i+1)+U_n1(j,i-1))+(h^2/4)*F(j,i)-(1/4)*u_y0(i+1);
        else
            U_n1(j,i) = (1/4)*(U_n1(j-1,i)+U_n1(j,i-1))+(h^2/4)*F(j,i)-(1/4)*(u_xb(j)+u_y0(i+1));
        end
    end
end

