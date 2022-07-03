%2D Heat Equation.
clear

%grid has n-2 interior points per dimension (overlapping)
nx = 20;  
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

% for stabiltiy fourier number should be less than 1/4  
r = D*dt/(dx.^2) ; % Fourier Number

TOL = 1e-6;

T = zeros(nx , nx , nt);
T(1,1:nx,:)  = 200;  %TOP
T(nx,1:nx,:) = 200;  %BOTTOM
T(1:nx,1,:)  = 200;  %LEFT
T(1:nx,nx,:) = 200;  %RIGHT

error = 1;
k = 0;

while error > TOL 
      k = k+1;
      Told = T;
      
      for i = 2:nx-1

          for j = 2:nx-1
              
              for l = 1 : nt                

              T(i,j,l) = (r)*((Told(i+1,j,l)-2*Told(i,j,l)+Told(i-1,j,l)) ... 
                      + (Told(i,j+1,l)-2*Told(i,j,l)+Told(i,j-1,l))) ...
                      + Told(i,j,l);                 
              end
          end

      end

      error = max(max(abs(Told-T)));
end

figure
for h = 1 : nt
    surf(x,y,T(:,:,h))
    hold on
end

 