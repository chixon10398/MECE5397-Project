clear all
close all
clc

%boundary lengths
a = 0;
b = 2*pi;
%discretize the problem
n = 100; %number of spaces
h = b/n; %distance between each position
x = 0:h:b;
y = x;
%create the functions
fb = y.*(b-y);
gb = ((b-y).^2).*cos(pi*y/b);
z = 1;
F = zeros(1,(n-1)^2); %create vector for F values
for j = h:h:b-h
    for i = h:h:b-h
        F(z) = sin(pi*i/b)*cos((pi/2)*(2*j/b+1));
        z = z+1;
    end
end
%Boundary conditions
u_x0 = fb;
u_xb = gb;
u_y0 = fb(1)+(x/b).*(gb(1)-fb(1));
%Build the vector for U values
U = zeros((n-1)^2,1);
U_D = ones(n-1,1); %assign a value for the Dirchielt boundary conditions

%ghost node
z=1;
for i = h:h:b-h
    F_D(z) = sin(pi*i/b)*cos(pi/2);
    z=z+1;
end

