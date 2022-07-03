function [ x, t, U ] = Crank_Nicolson( vString, fString, a, N, M,g1,g2 )
% The Crank Nicolson provides a solution to the parabolic equation provided
%   The Crank Nicolson method uses linear system of equations to solve the
%   parabolic equation.

%Prepare the grid and grid spacing variables. 
dt = 1/M;
t = dt * [0:M];
h = 1/N;
x = 2 + h * [0:N]';%Shift x by 2 that way we have 2 <= x <= 3

%Prepare the matrix that will store the solutions over the grid
U = zeros(N+1, M+1);
%This will fill the first column with the initial condition. feval will
%evaluate the initial condition at all values of x
U(:,1) = feval(vString, x);
%This fills the boundary conditions. One boundary condition goes on the
%first row the other boundary condition goes on the last row
U(1,:) = feval(g1, t);
U(end,:) = feval(g2, t);

%The loop that will populate the matrix with the solution
n = 1:N+1;%Start at 2 since n=1 is the initial condition
e = ones(N+1,1);
B = spdiags([-1*e 2*e -1*e],-1:1, N+1, N+1)*(1/h^2);
A = (speye(N+1)+((a*dt)/2)*B);
X = (speye(N+1)-((a*dt)/2)*B);
R = chol(A);%Choleski decomposition
for m=2:M+1
    %The linear system is solved. 
    b = X*U(n,m-1) + dt * feval(fString, x(n), (t(m)+t(m-1))*0.5);
    b = R'\b;
    U(n,m) = R\b;
end
end