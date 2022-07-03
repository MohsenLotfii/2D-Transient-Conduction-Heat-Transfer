%=========================================================================%
%========== IMPLICIT CRANK-NICOLSON METHOD FOR 2D HEAT EQUATION ==========%
% . . .  Implicit solution of the 2D Heat Equation : Crank-Nicolson . . .
%  Solves the 2D heat equation with an implicit finite difference scheme
%  Created by: 'Mohsen Lotfi'
%  Last Revised: 10 May 2017
%  Persian Gulf University, Bushehr, Iran
%=========================================================================%

% tic % Start stopwatch timer

% clear workspace
clear

% --------------------------- Physical parameters ------------------------%
nx = 5;
ny = 5;
x = linspace(0,1,nx);
dx = x(2)-x(1);
y = x;
dy = dx;

%time interval
tf = 10;
nt = 20;
dt = tf/nt;

% Properties of Copper
k_s = 387             ; %W/(m K)  Thermal conductivity - conduction coefficient
rho = 8940            ; %kg.m^3   Density
Cp = 380              ; %J/kg K   Spefic heat constant pressure
D = k_s/(rho*Cp)      ; %m^2/s    Diffusivity

r = D*dt/(dx.^2) ;

% ---------------- constructing penta-diagonal matrix P ------------------%
% for diagonal elements
for i = 1 : nx*ny
    P(i,i)= 2*((1/dx^2)+(1/dy^2));
end

% for uper & lower diagonal elements
for i = 1 : nx-1
    for j = 1 : ny
        P(i+(j-1)*nx , i+(j-1)*nx+1)= -1/dx^2;
        P(i+(j-1)*nx+1 , i+(j-1)*nx)= -1/dx^2;
    end
end

% for off-diagonal elements
for i = 1 : nx
    for j = 1 : ny-1
        P(i+(j-1)*nx , i+j*nx)= -1/dy^2;
        P(i+j*nx , i+(j-1)*nx)= -1/dy^2;
    end
end

I = eye(nx*ny , nx*ny) ;
H = I - (r/2)*P ;

% ?????? ??? ???? ???%

% ------------------- solving the equation A*x = B -----------------------%

% Set Initial Condition & Boundary Conditions
T = 100 + zeros(nx , ny , nt); % initial temperature Ti=100;
T(1,1:nx,:)  = 200;  %TOP
T(nx,1:nx,:) = 200;  %BOTTOM
T(1:nx,1,:)  = 200;  %LEFT
T(1:nx,nx,:) = 200;  %RIGHT

% This is the time advance loop.
Time = 0;
for n = 1 : nt
    
    % we can put the knowns in a vector B
    B = zeros(nx*ny , nx*ny , 1);
    
    % Compute B vector:
    for i = 1:ny
        for j = 1:nx
            
            B(2:nx*ny-1 , 2:nx*ny-1) = (I-(r/2)*P)* Told(:,:,n);
            B(1 , 1) = 200;
            B(nx*ny , nx*ny) = 200;

        end
    end

    % Compute solution vector (list of the unknowns)
    Tnew_vector = A\B;
    
    % Create 2D matrix from vector
    Tnew = Tnew_vector(n);
    T = Tnew;
    Time = Time + dt;
    
end

% toc % Read elapsed time from stopwatch 

