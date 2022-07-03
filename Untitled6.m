
a = -2*pi;
b =  2*pi;

NN = 2.^(3:9);
for jj = 1 : length(NN)
    n = NN(jj);               %Number of grid points
    dx = (b-a)/(n-1);         %Spatial mesh spacing
    x = a : dx : b;           %Mesh in x direction
    c = 1;
    u = @(x,t) cos(c*t).*cos(x);
    u0 = u(x,0.0)';

    e = ones(n,1);
    A = spdiags([e -2*e e], -1:1, n, n);   %Laplacian matrix
    A(1,end) = 1;  A(end,1) = 1;           %Periodic boundary conditions

    T = 1;                                 %Final simulation time
    dt = 1e-3;                             %Temporal mesh spacing
    nsteps = (floor(T/dt));                %Number if time steps


    %u_old = u0;  %First order approxiation of u_t(x,0)=0.
    %u_now = u0; 

                  %Second order approxiation of u_t(x,0)=0.
    u_old = u0;
    u_now = u0 + 0.5*(dt*c/dx)^2*(A* u0);  %(approximates u(x,dt))

    for ii = 2 : nsteps  %First time step was directly above
       u_new = 2*u_now - u_old + (dt*c/dx)^2*(A* u_now);
       u_old = u_now;
       u_now = u_new;
    end

    LI(jj) =  norm(u_new - u(x,T)',inf);
    L2(jj) = dx*norm(u_new - u(x,T)',2);
    if ( jj > 1 )  %Look at order of convergence
       log ( L2(jj-1)/L2(jj) ) / log ( 2 ) 
    end
end

figure(1)
semilogy( L2,'-')
figure(2)
plot(x,u_new,'ro-',x,u(x,T),'--')

