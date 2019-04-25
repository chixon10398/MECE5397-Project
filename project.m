clear all
close all
clc

%boundary lengths
a = 0;
b = 2*pi;
%discretize the problem
n = 100; %number of spaces
h = b/n; %distance between each position
x = h:h:b;
y = x;
%create the functions
fb = y.*(b-y);
gb = ((b-y).^2).*cos(pi*y/b);
F = sin(pi*x/b).*cos((pi/2)*(2*y/b+1));
%Boundary conditions
u_x0 = fb;
u_xb = gb;
u_y0 = fb(1)+(x/b).*(gb(1)-fb(1));
%ghost node
