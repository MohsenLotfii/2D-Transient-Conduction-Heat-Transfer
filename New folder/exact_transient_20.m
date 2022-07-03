
tic % Start stopwatch timer
    
% clear work space & command window
clear, clc

% Boundary Conditions & Initial Condition
T0 = 200; 
Ti = 100;

L = 1; % Width of domain in x direction
H = 1; % Height of domain in y direction

k_s = 387;        %W/(m K)  Thermal conductivity - conduction coefficient
rho = 8940;       %kg.m^3   Density
Cp = 380;         %J/kg K   Spefic heat constant pressure
D = k_s/(rho*Cp); %m^2/s    Diffusivity

a = 0; 
b = 1;
N = 20;
h = (b-a)/(N-1); % Mesh spacing

t = 50;
% del_t = 0.01;
% t_range = 0 : del_t : t ;

x_range = 0:h:1;
y_range = 0:h:1;

% x = [0:N-1];
% y = [0:N-1];
% X = repmat(x,N,1);
% Y = repmat(y(:),1,N);

for i = 1:length(x_range)
    x = x_range(i);

    for j = 1:length(y_range)
        y = y_range(j);

        Zigma = 0;
        
        for n = 1 : 60
            beta = n*pi/L;
            
            for m = 1 : 60
                gama = m*pi/H;
                F_mn = (4*(Ti-T0)/(L*H))*(1/gama)*(1/beta)*((-1)^n-1)*((-1)^m-1);
                Zigma = Zigma + F_mn*((sin(beta*x))*exp(-D*t*((beta)^2))*(sin(gama*y))*exp(-D*t*((gama)^2)));
            end
            
        end
        
        T(i,j)= Zigma + T0 ;

    end % for j = 1:length(y_range)

end % for i = 1:length(x_range)

T_exact = T;
 
% Plotting the exact solution
figure(1)
%[X,Y] = meshgrid(x_range,y_range)
%contourf(X,Y,T_exact)
contourf(x_range,y_range,T_exact)
title('Temperature in t = 121.5830s')
xlabel('x')
ylabel('y')
colorbar

save T_exact

toc % Read elapsed time from stopwatch

