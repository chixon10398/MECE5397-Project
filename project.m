clear all
close all
clc
%initial run
%boundary lengths
a = 0;
b = 2*pi;
%discretize the problem
n = 25; %number of steps to b
h = b/n; %distance between each position
x = 0:h:b;
y = b:-h:0;
%tolerance for acceptable value
tol = 10^-5;
%create the functions
fb = y.*(b-y).^2;
gb = ((b-y).^2).*cos(pi*y/b);
F = zeros(n-1,n-1); %create vector for F values
%% new test
GaussS = zeros(n-1,n-1); %U for n
GaussSn1 = zeros(n-1,n-1); % U for n+1
SOR = zeros(n-1,n-1);
SORn1 = zeros(n-1,n-1);
save('variables','GaussS','GaussSn1')
save('variables','SOR','SORn1')
%% Gauss Seidel Method Checkpoint
load('variables','GaussS','GaussSn1')

ghost = zeros(1,n-1);
g = ghost;
%create the elements for F
for i = 1:n-1
    for j = 1:n-1
        F(j,i) = sin(pi*((x(i+1)-a)/(b-a)))*cos((pi/2)*((2*((y(j+1)-a)/(b-a)))+1));
    end
end

%Boundary conditions
u_x0 = fb;
u_xb = gb;
u_y0 = fb(n+1)+(x/b).*(gb(n+1)-fb(n+1));

index = 1;
diff = 1;
while diff > tol
    for i = 1:n-1 %this loop adds the ghost boundarys
        if i == 1
            ghost(i) = (1/4)*(2*GaussS(1,i)+g(i+1)+u_x0(n+1));
            GaussSn1(1,i) = (1/4)*(GaussS(2,i)+GaussS(1,i+1)+h^2*F(1,i)+ghost(i)+u_x0(2));
        elseif i == n-1
            ghost(i) = (1/4)*(2*GaussS(1,i)+g(i-1)+u_xb(n+1));
            GaussSn1(1,i) = (1/4)*(GaussS(2,i)+GaussSn1(j,i-1)+h^2*F(j,i)+ghost(i)+u_xb(2));
        else
            ghost(i) = (1/4)*(2*GaussS(1,i)+g(i+1)+g(i-1));
            GaussSn1(1,i) = (1/4)*(GaussS(2,i)+GaussS(j,i+1)+GaussSn1(j,i-1)+h^2*F(j,i)+ghost(i));
        end
    end
    g = ghost;
    for j = 2:n-1
        for i = 2:n-1
            if j == n-1 && i == n-1
                GaussSn1(j,i) = (1/4)*(GaussSn1(j-1,i)+GaussSn1(j,i-1)+h^2*F(j,i)+u_xb(j+1)+u_y0(i+1));
            elseif j == n-1 && i < n-1
                GaussSn1(j,i) = (1/4)*(GaussSn1(j-1,i)+GaussS(j,i+1)+GaussSn1(j,i-1)+h^2*F(j,i)+u_y0(i+1));
            elseif j < n-1 && i == n-1
                GaussSn1(j,i) = (1/4)*(GaussS(j+1,i)+GaussSn1(j-1,i)+GaussSn1(j,i-1)+h^2*F(j,i)+u_xb(j+1));
            elseif j < n-1 && i < n-1
                GaussSn1(j,i) = (1/4)*(GaussS(j+1,i)+GaussSn1(j-1,i)+GaussS(j,i+1)+GaussSn1(j,i-1)+h^2*F(j,i));
            end
        end
        if j == n-1
            GaussSn1(j,1) = (1/4)*(GaussSn1(j-1,1)+GaussS(j,2)+h^2*F(j,1)+u_x0(j+1)+u_y0(2));
        elseif j < n-1
            GaussSn1(j,1) = (1/4)*(GaussS(j+1,1)+GaussSn1(j-1,1)+GaussS(j,2)+h^2*F(j,i)+u_x0(j+1));
        end
    end
save('variables','GaussSn1')
diff = norm(GaussSn1-GaussS);
GaussS = GaussSn1;
index = index+1
save('variables','GaussS')
end
figure(1)
surf(GaussS)

%% SOR method check point
load('variables','SOR','SORn1')
w = 1.1;
%create the elements for F
for i = 1:n-1
    for j = 1:n-1
        F(j,i) = sin(pi*((x(i+1)-a)/(b-a)))*cos((pi/2)*((2*((y(j+1)-a)/(b-a)))+1));
    end
end

%Boundary conditions
u_x0 = fb;
u_xb = gb;
u_y0 = fb(n+1)+(x/b).*(gb(n+1)-fb(n+1));

index2 = 1;
diff = 1;
while diff > tol
    for i = 1:n-1 %this loop adds the ghost boundarys
        if i == 1
            SORn1(1,i) = (1-w)*SOR(1,i)+(w/4)*(SOR(2,i)+SOR(1,i+1)+h^2*F(1,i)+ghost(i)+u_x0(2));
        elseif i == n-1
            SORn1(1,i) = (1-w)*SOR(1,i)+(w/4)*(SOR(2,i)+SORn1(j,i-1)+h^2*F(j,i)+ghost(i)+u_xb(2));
        else
            SORn1(1,i) = (1-w)*SOR(1,i)+(w/4)*(SOR(2,i)+SOR(j,i+1)+SORn1(j,i-1)+h^2*F(j,i)+ghost(i));
        end
    end
    for j = 2:n-1
        for i = 2:n-1
            if j == n-1 && i == n-1
                SORn1(j,i) = (1-w)*SOR(j,i)+(w/4)*(SORn1(j-1,i)+SORn1(j,i-1)+h^2*F(j,i)+u_xb(j+1)+u_y0(i+1));
            elseif j == n-1 && i < n-1
                SORn1(j,i) = (1-w)*SOR(j,i)+(w/4)*(SORn1(j-1,i)+SOR(j,i+1)+SORn1(j,i-1)+h^2*F(j,i)+u_y0(i+1));
            elseif j < n-1 && i == n-1
                SORn1(j,i) = (1-w)*SOR(j,i)+(w/4)*(SOR(j+1,i)+SORn1(j-1,i)+SORn1(j,i-1)+h^2*F(j,i)+u_xb(j+1));
            elseif j < n-1 && i < n-1
                SORn1(j,i) = (1-w)*SOR(j,i)+(w/4)*(SOR(j+1,i)+SORn1(j-1,i)+SOR(j,i+1)+SORn1(j,i-1)+h^2*F(j,i));
            end
        end
        if j == n-1
            SORn1(j,1) = (1-w)*SOR(j,i)+(w/4)*(SORn1(j-1,1)+SOR(j,2)+h^2*F(j,1)+u_x0(j+1)+u_y0(2));
        elseif j < n-1
            SORn1(j,1) = (1-w)*SOR(j,i)+(w/4)*(SOR(j+1,1)+SORn1(j-1,1)+SOR(j,2)+h^2*F(j,i)+u_x0(j+1));
        end
    end
save('variables','SORn1')
diff = norm(SORn1-SOR);
SOR = SORn1;
index2 = index2+1
save('variables','SOR')
end

figure(2)
surf(SOR)
