%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by : Mohsen Lotfi
% last modified : 5 May 2017

% Clear memory, close figures, clear command window for every run
clear
clc

 syms r T 
 % Define input variables
L_i = 1               ; %m  length in verticle direction
L_j = 1               ; %m  length in horizontal direction
n_i = 20              ; %   Number of nodes in i direction
n_j = 20              ; %   Number of nodes in j direction (adjusted to keep constant spacing)
del_x = L_i/(n_i-1)   ; %m  Spacing of nodes

nt = 500; % Number of timesteps to compute

% Properties of Copper
k_s = 387             ; %W/(m K)  Thermal conductivity - conduction coefficient
rho = 8940            ; %kg.m^3   Density
Cp = 380              ; %J/kg K   Spefic heat constant pressure
D = k_s/(rho*Cp)  ; %m^2/s    Diffusivity

% Compute stable timestep
% del_t = del_x^2/D/4;
del_t = 0.1;

% for stabiltiy fourier number should be less than 1/4  
r = D*del_t/(del_x.^2) ; % Fourier Number

k = 0 ; % iteration counter
n = 100;
% define T array
T = zeros(n);
T(1,1:n) = 200;  %TOP
T(n,1:n) = 200;  %BOTTOM
T(1:n,1) = 200;  %LEFT
T(1:n,n) = 200;  %RIGHT

Tnew = T; % T_k+1

%  time = 0;
 
 
for n=0:del_t:nt
    
    % update iteration counter 
%     k = k + 1;
    
    for i = 2:del_x:n_i
        
        for j = 2:del_x:n_j
            
            % index: n+1
            Tnew(i,j,n+1) = T(i,j,n)+r*((T(i,j+1,n)-2*T(i,j,n)+T(i,j-1,n))+(T(i+1,j,n)-2*T(i,j,n)+T(i-1,j,n)));
            
        end % (j loop)
        
    end % (i loop)

    % update T
    T = Tnew ;
end

Temp = T;
time = time+del_t;

Temp(1:n_i,1:n_j)= T(2:n_i+1,1:n_j,k-1) ;

%Creating x direction position matrix
for i = 1:n_i
    x(i) = del_x*(i-1) ;
end

%Creating y direction position matrix
for j = 1:n_j
    y(j)= del_x*(j-1) ;
end

figure(13)
surf(y,x,Temp(:,:))      % surf(column matrix, row matrix, 2D temp matrix)

