%=================== 2D Heat Equation (Parabolic Problem) ================%
clear all
close all
clc

% grid has n-2 interior points per dimension (overlapping)
nx = 20;  
x = linspace(0,1,nx);
dx = x(2)-x(1);
y = x;
dy = dx;

% time interval
tf = 10; % end time
nt = 20; % number of time steps
dt = tf/nt;

% Properties of Copper
k_s = 387         ; %W/(m K)  Thermal conductivity - conduction coefficient
rho = 8940        ; %kg.m^3   Density
Cp = 380          ; %J/kg K   Spefic heat constant pressure
D = k_s/(rho*Cp)  ; %m^2/s    Diffusivity

% for stabiltiy fourier number should be less than 1/4  
r = (D*dt)/(dx^2) ; % Fourier Number

% Set IC & BC
Ti = 100; % initial temperature Ti=100;
T = Ti + zeros(nx , nx); 
T(1,1:nx)  = 200;  %TOP
T(nx,1:nx) = 200;  %BOTTOM
T(1:nx,1)  = 200;  %LEFT
T(1:nx,nx) = 200;  %RIGHT

TOL = 1e-6;
error = 1;

k = 0; %inital iteration step number


while error > TOL %m aximum allowable error
      k = k+1; 
      Told = T;
      
      for i = 2:nx-1

          for j = 2:nx-1         

              T(i,j) =  Told(i,j) ...
                  + r*(Told(i-1,j)-2*Told(i,j)+Told(i+1,j) ...
                  + Told(i,j-1)-2*Told(i,j)+Told(i,j+1));          

          end % for j = 2:nx-1

      end % for i = 2:nx-1

      error = max(max(abs(Told-T)));
end

figure(2)
contourf(x,y,T)
title('Temperature in t=10s')
xlabel('x')
ylabel('y')
colorbar

