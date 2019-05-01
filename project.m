clear all
close all
clc

%boundary lengths
a = 0;
b = 2*pi;
%discretize the problem
n = 10; %number of steps to b
h = b/n; %distance between each position
x = 0:h:b;
y = x;
%tolerance for acceptable value
tol = 10^-5;
%create the functions
fb = y.*(b-y).^2;
gb = ((b-y).^2).*cos(pi*y/b);
F = zeros(n-1,n-1); %create vector for F values
U_n = zeros(n-1,n-1); %U for n
U_n1 = zeros(n-1,n-1); % U for n+1
save('variables','U_n','U_n1')
load('variables','U_n','U_n1')
ghost = zeros(1,n-1); %wher the values for our ghost node will be stored
%create the elements for F
for i = 1:n-1
    for j = 1:n-1
        F(j,i) = sin(pi*((x(i+1)-a)/(b-a)))*cos((pi/2)*((2*((y(j+1)-a)/(b-a)))+1));
    end
end

%Boundary conditions
u_x0 = fb;
u_xb = gb;
u_y0 = fb(1)+(x/b).*(gb(1)-fb(1));

%ghost node
for i = 1:n-1
    if i == 1
        ghost(i) = (1/4)*(2*U_n(i,2)+U_n(i+1,1)+h^2*F(i,1)+u_x0(i));
    elseif i == n-1
        ghost(i) = (1/4)*(2*U_n(i,2)+U_n(i-1,1)+h^2*F(i,j)+u_xb(i));
    else
        ghost(i) = (1/4)*(2*U_n(i,2)+U_n(i+1,1)+U_n(i-1,1)+h^2*F(j,i));
    end
end

%Iterate for the next values
index = 1;
diff = 1;
while diff > tol
for j = 1:n-1
    for i = 1:n-1
        if j == 1 && i == 1
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n(j,i+1)+h^2*F(j,i)+ghost(i)+u_x0(j+1));
        elseif j == 1 && i == n-1
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n1(j,i-1)+h^2*F(j,i)+ghost(i) + u_xb(j+1));
        elseif j == n-1 && i == 1
            U_n1(j,i) = (1/4)*(U_n1(j-1,i)+U_n(j,i+1)+h^2*F(j,i)+u_x0(j+1)+u_y0(i+1));
        elseif j == n-1 && i == n-1
            U_n1(j,i) = (1/4)*(U_n1(j-1,i)+U_n1(j,i-1)+h^2*F(j,i)+u_xb(j+1)+u_y0(i+1));
        elseif j == 1 && i > 1 && i < n-1
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n(j,i+1)+U_n1(j,i-1)+h^2*F(j,i)+ghost(i));
        elseif j == n-1 && i > 1 && i < n-1
            U_n1(j,i) = (1/4)*(U_n1(j-1,i)+U_n(j,i+1)+U_n1(j,i-1)+h^2*F(j,i)+u_y0(i+1));
        elseif i == 1 && j > 1 && j < n-1
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n1(j-1,i)+U_n(j,i+1)+h^2*F(j,i)+u_x0(j+1));
        elseif i == n-1 && j > 1 && j < n-1
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n1(j-1,i)+U_n1(j,i-1)+h^2*F(j,i)+u_xb(j+1));
        else            
            U_n1(j,i) = (1/4)*(U_n(j+1,i)+U_n1(j-1,i)+U_n(j,i+1)+U_n1(j,i-1)+h^2*F(j,i));
        end 
        save('variables','U_n1')
    end
end
diff = norm(U_n1-U_n);
U_n = U_n1;
index = index+1
save('variables','U_n')
end
