% %=================== 2D Heat Equation (Parabolic Problem) ================%
% clear
% 
% % grid has n-2 interior points per dimension (overlapping)
% nx = 20;  
% x = linspace(0,1,nx);
% dx = x(2)-x(1);
% y = x;
% dy = dx;
% 
% % time interval
% tf = 10; % end time
% nt = 20; % number of time steps
% dt = tf/nt;
% 
% % Properties of Copper
% k_s = 387         ; %W/(m K)  Thermal conductivity - conduction coefficient
% rho = 8940        ; %kg.m^3   Density
% Cp = 380          ; %J/kg K   Spefic heat constant pressure
% D = k_s/(rho*Cp)  ; %m^2/s    Diffusivity
% 
% % for stabiltiy fourier number should be less than 1/4  
% r = (D*dt)/(dx^2) ; % Fourier Number
% 
% % Set IC & BC
% Ti = 100; % initial temperature Ti=100;
% T = Ti + zeros(nx , nx , nt); 
% T(1,1:nx,:)  = 200;  %TOP
% T(nx,1:nx,:) = 200;  %BOTTOM
% T(1:nx,1,:)  = 200;  %LEFT
% T(1:nx,nx,:) = 200;  %RIGHT
% 
% TOL = 1e-6;
% error = 1;
% k = 0; %inital iteration step number
% 
% while error > TOL %m aximum allowable error
%       k = k+1;
%       Told = T;
%       
%       if dt <= (dx^2)/(2*D)
%       
%           for i = 2:nx-1
% 
%               for j = 2:nx-1
% 
%                   for l = 1 : nt % Loop over time steps           
% 
%                   T(i,j,l) =  Told(i,j,l)...
%                       + r*(Told(i-1,j,l)-2*Told(i,j,l)+Told(i+1,j,l)...
%                       + Told(i,j-1,l)-2*Told(i,j,l)+Told(i,j+1,l));
%                
%                   end % for l = 1 : nt
% 
%               end % for j = 2:nx-1
% 
%           end % for i = 2:nx-1
%       
%       end
% 
%       error = max(max(abs(Told-T)));
% end
% 
% figure(1)
% for h = 1 : nt
%     surf(x,y,T(:,:,h))
%     hold on
% end
% title('2D Transient Heat Transfer')
% xlabel('x')
% ylabel('y')
% zlabel('Temperature')
% colorbar
% 
% figure(2)
% contourf(x,y,T(:,:,nt))
% title('Temperature in t=10s')
% xlabel('x')
% ylabel('y')
% colorbar
% 

%=================== 2D Heat Equation (Parabolic Problem) ================%
clear

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
% r = (D*dt)/(dx^2) ; % Fourier Number

% Set IC & BC
Ti = 100; % initial temperature Ti=100;
% T = Ti + zeros(nx , nx , nt); 
% T(1,1:nx,:)  = 200;  %TOP
% T(nx,1:nx,:) = 200;  %BOTTOM
% T(1:nx,1,:)  = 200;  %LEFT
% T(1:nx,nx,:) = 200;  %RIGHT

T = 0;
TOL = 1e-6;
error = 1;
k = 0; %inital iteration step number
t = 0;

while error > TOL %m aximum allowable error
      k = k+1;
      Told = T;
      
      if dt <= (dx^2)/(2*D)
          
       for l = 1 : nt % Loop over time steps
           
                 T = zeros(nx,nx);
                 r = (D*dt)/(dx^2) ; % Fourier Number
     
          for i = 2:nx-1

                for j = 2:nx-1                 
                 
                  T(i,j) =  Told(i,j)...
                      + r*(Told(i-1,j)-2*Told(i,j)+Told(i+1,j)...
                      + Told(i,j-1)-2*Told(i,j)+Told(i,j+1));
               
                end % for l = 1 : nt

           end % for j = 2:nx-1
              
       end % for i = 2:nx-1
          
                T(1,1:nx)  = 200;  %TOP
                T(nx,1:nx) = 200;  %BOTTOM
                T(1:nx,1)  = 200;  %LEFT
                T(1:nx,nx) = 200;  %RIGHT
                
       
      end

      
            t = t+dt;
%       error = max(max(abs(Told-T)));
end

figure(1)
for h = 1 : nt
    surf(x,y,T(:,:,h))
    hold on
end
title('2D Transient Heat Transfer')
xlabel('x')
ylabel('y')
zlabel('Temperature')
colorbar

figure(2)
contourf(x,y,T(:,:,nt))
title('Temperature in t=10s')
xlabel('x')
ylabel('y')
colorbar

