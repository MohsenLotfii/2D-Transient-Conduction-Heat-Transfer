%=========================================================================%
%========== IMPLICIT CRANK-NICOLSON METHOD FOR 2D HEAT EQUATION ==========%
% . . .  Implicit solution of the 2D Heat Equation : Crank-Nicolson . . .
%  Solves the 2D heat equation with an implicit finite difference scheme
%  Created by: 'Mohsen Lotfi'
%  Last Revised: 10 May 2017
%  Persian Gulf University, Bushehr, Iran
%=========================================================================%

tic % Start stopwatch timer
% clear workspace
clear

% --------------------------- Physical parameters ----------------------- %
% time interval
tf = 2000;              % end time
nt = 400;              % number of time steps
dt = tf/nt;           % time spacing

% space grids
L = 1;                % length and width of plate
nx = 20;              % number of elements in x-direction
ny = 20;              % number of elements in x-direction
x = linspace(0,1,nx); % vector of x values
y = x;                % vector of y values
dx = x(2)-x(1);       % grid spacing
h = dx;

[xx , yy] = meshgrid(x,y);

% Properties of Copper
k_s = 387             ; %W/(m K)  Thermal conductivity - conduction coefficient
rho = 8940            ; %kg.m^3   Density
Cp = 380              ; %J/kg K   Spefic heat constant pressure
k = k_s/(rho*Cp)      ; %m^2/s    Diffusivity 

% ---------------- constructing penta-diagonal matrix H ------------------%
% set constants a and c
a = (1 + (2*k)/(h^2));
c = (-k/(2*(h^2)));
d = (1 - (2*k)/(h^2));
e = (k/(2*(h^2)));

% set the matrix for H
A = a*eye(nx) + c*diag(ones(1, nx-1), 1) + c*diag(ones(1, nx-1),-1);
B = c*eye(nx);
D = d*eye(nx) + e*diag(ones(1, nx-1), 1) + e*diag(ones(1, nx-1),-1);
E = e*eye(nx);
C = zeros(nx);

% Create the matrix H, composed of sub-matrices A,B and C
Ap = eye(nx);
Bp = diag(ones(1, nx-1), 1) + diag(ones(1, nx-1), -1);
Dp = eye(nx);
Ep = diag(ones(1, nx-1), 1) + diag(ones(1, nx-1), -1);
Cp = ones(nx) - Ap - Bp;

% Define the H_1 matrix using KRON to "replace" 1's in Ap, Bp, and Cp
%with copies of A, B, and C
H_1 = kron(Ap, A) + kron(Bp, B) + kron(Cp, C);

% Define the H_2 matrix using KRON to "replace" 1's in Ap, Bp, and Cp
%with copies of A, B, and C
H_2 = kron(Dp, D) + kron(Ep, E) + kron(Cp, C);

% ----------- Set the Initial Condition & Boundary Conditions ------------%

% I.C.
i_c = 100; % the initial cond value for the plate
T = 100 * ones(nx,nx); % initial temperature Ti=100;

% B.C.
top = 200;
T(1,:) = top;          % TOP
left = 200;
T(:,1) = left;         % LEFT
bottom = 200;
T(nx,:) = bottom;      % BOTTOM
right = 200;
T(:,nx) = right;       % RIGHT

% ---------------------------- Eq. A*x = B -------------------------------%
% There are several standard algorithms to solve this equation,
% given in the following table:

% Type = D means the algorithm is direct,...
% i.e. it gives the exact answer after the indicated number of steps.

% Type = I means the algorithm is iterative,
% i.e. that it repeatedly updates an approximate solution,...
% so that the solution converges as we do more iterations.
%-------------------------------------------------------------------------%
%    Algorithm   Type  Serial Time     PRAM Time      Storage  #Procs 
%    ---------   ----  -----------     ---------      -------  ------
%    Dense LU      D       N^3              N           N^2      N^2
%    Band LU       D       N^2              N         N^(3/2)     N
%    Jacobi        I       N^2              N            N        N
%    Inv(H)*b      D       N^2            log N         N^2      N^2
%    CG            I       N^(3/2)      N^(1/2)*log N    N        N
%    SOR           I       N^(3/2)      N^(1/2)          N        N
%    Sparse LU     D       N^(3/2)      N^(1/2)       N*log N     N
%    FFT           D       N*log N        log N          N        N
%    Multigrid     I       N              (log N)^2      N        N
%    Lower Bound           N              log N          N 
% 
% Key to abbreviations:
%      Dense LU: Gaussian elimination, treating H as dense
%      Band LU : Gaussian elimination, treating H as zero outside a band
%                                      of half-width n-1 near diagonal
%      Sparse LU : Gaussian elimination, exploiting entire 
%                                        zero-structure of H
%      Inv(H)*b : precompute and store inverse of H, multiply it by 
%                   right-hand-side b
%      CG      : Conjugate Gradient method
%      SOR     : Successive Overrelaxation
%      FFT     : Fast Fourier Transform based method
%-------------------------------------------------------------------------%

T = T(:);

[Low,Up,Per] = lu(H_1); % LU matrix factorization

% The solution to Ax = b is obtained with matrix division x = A\b
% The solution is actually computed by solving two triangular systems:
% y = L\b  =======>  x = U\y

% This is the time advance loop.
for m = 1 : dt : nt

    % we can put the knowns in a vector B
    % Compute B vector:
    B = H_2*T;
    
    % Compute solution vector (list of the unknowns)
    T = Up\(Low\(Per*B));
    
    % replace the B.C. values
    T(1:nx) = left;                    % Left B.C.
    T(((2*nx):nx:nx*nx)) = bottom;     % Bottom B.C.
    T((nx+1):nx:(nx*nx)) = top;        % Top B.C.
    T((((nx-1)*nx)+1):(nx^2)) = right; % Right B.C.    
end

% create a surf plot of the temperature (y-axis) against grid positions
TT = reshape(T,nx,nx);

figure(4)
surf(xx,yy,TT)
title('2D Transient Heat Transfer in t = 800 sec. and \Deltat = 5 sec.')
xlabel('x')
ylabel('y')
zlabel('Temperature')
colorbar

figure(5)
contourf(xx,yy,TT(:,:))
title('Temperature in in t = 800 sec. and \Deltat = 5 sec.')
xlabel('x')
ylabel('y')
colorbar

% Calculating Error:
load T_exact

Error = abs(T_exact - T);
Err = (norm(Error)/norm(T_exact))

toc % Read elapsed time from stopwatch 

