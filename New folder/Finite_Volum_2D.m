% Create Grid and number of cells
a = 0; 
b = 1; 
N = 10;
% Define edges
x_edges = linspace(a,b,N+1);
y_edges = linspace(a,b,N+1);
%Define distance between edges
delta_x = x_edges(2) - x_edges(1);
delta_y = y_edges(2) - y_edges(1);
%Define cell centers
x_centers = a+delta_x/2 : delta_x : b;
y_centers = a+delta_y/2 : delta_y : b;
[X,Y] = meshgrid(x_centers,y_centers);      % 2d arrays of x,y center values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % (i,j) cell.
%Initialize solution array
U = zeros(N,N);
 
% Fill the solution array with initial data
U = X*0+1;
%define distance function
phi = sqrt((X-1/2).^2+(Y-1/2).^2);
gradphix = (X-1/2)./phi;
gradphiy = (Y-1/2)./phi;
% %reverse sign so that the gradient points inwards
v1 = -gradphix;
v2 = -gradphiy;
 
%Define timestep
delta_t = 0.001;
F = v1.*U;
G = v2.*U;
%Store initial data for later use
ic = U;
timesteps = 20;
for k = 1:timesteps
for j = 1:N % denote row position
 
     %denote column position
     for i = 1:N
 
        %Top row
        if (j ==N) 
            %Down flux
            G_down =  (1/delta_y)*(G(i,j)-G(i,j-1));
 
            %Nothing is coming down from above
            G_up =  0;
 
        %Bottom row    
        elseif (j == 1)
 
            %Nothing is coming in from below
            G_down = 0;
 
            %Up flux
            G_up = (1/delta_y)*(G(i,j+1)-G(i,j));
 
        %Rows in the middle
        else
 
            G_up = (1/delta_y)*(G(i,j+1)-G(i,j));
            G_down =  (1/delta_y)*(G(i,j)-G(i,j-1));
 
        end
 
 
        %Sideway fluxes
            if i == 1 
                F_left = 0;
                F_right = (1/delta_x)*(F(i+1,j)-F(i,j));
            elseif i == N
                F_left = (1/delta_x)*(F(i,j)-F(i-1,j));
                F_right = 0;
            else
                F_right = (1/delta_x)*(F(i+1,j)-F(i,j));
                F_left = (1/delta_x)*(F(i,j)-F(i-1,j));
            end
 
        Unew(i,j) = U(i,j) - (delta_t)*(F_right-F_left)-(delta_t)*(G_up - ...
             G_down);
     end
end
U = Unew;
surf(X,Y,Unew)
pause(0.1);
end

% % From Michael
% % 1)Your main error is the flux ! see below and compare
% % 2) Secondly use reflection b.c.'s.
% % 3) The scheme below is actually FTCS which is unstable for convection
% % alone !!! so beware...
% % Let me know ! what you is your project ?
% % Create Grid and number of cells
% a = 0;
% b = 1;
% N = 40;
% % Define edges
% x_edges = linspace(a,b,N+1);
% y_edges = linspace(a,b,N+1);
% %Define distance between edges
% delta_x = x_edges(2) - x_edges(1);
% delta_y = y_edges(2) - y_edges(1);
% %Define cell centers
% x_centers = a+delta_x/2 : delta_x : b;
% y_centers = a+delta_y/2 : delta_y : b;
% [X,Y] = meshgrid(x_centers,y_centers); % 2d arrays of x,y center values
% X = X'; % transpose so that X(i,j),Y(i,j) are
% Y = Y'; % (i,j) cell.
% %Initialize solution array
% U = zeros(N,N); Unew = zeros(N,N);
% 
% % Fill the solution array with initial data
% U = X*0+1;
% %define distance function
% phi = sqrt((X-1/2).^2+(Y-1/2).^2);
% gradphix = (X-1/2)./phi;
% gradphiy = (Y-1/2)./phi;
% % %reverse sign so that the gradient points inwards
% v1 = -gradphix;
% v2 = -gradphiy;
% 
% %Define timestep
% delta_t = 0.001;
% F = v1.*U;
% G = v2.*U;
% %Store initial data for later use
% ic = U;
% timesteps = 20;
% for k = 1:timesteps
% for j = 1:N % denote row position
% 
% %denote column position
% for i = 1:N
% 
% %Top row
% if (j ==N)
% %Down flux
% G_down = (1/delta_y)*0.5*(G(i,j)+G(i,j-1));
% 
% %Nothing is coming down from above
% % G_up = 0;
% G_up = G_down;
% 
% %Bottom row
% elseif (j == 1)
% 
% %Nothing is coming in from below
% % G_down = 0;
% 
% %Up flux
% G_up = (1/delta_y)*0.5*(G(i,j+1)+G(i,j));
% G_down =G_up;
% %Rows in the middle
% else
% 
% G_up = (1/delta_y)*0.5*(G(i,j+1)+G(i,j));
% G_down = (1/delta_y)*0.5*(G(i,j)+G(i,j-1));
% 
% end
% 
% 
% %Sideway fluxes
% if i == 1
% % F_left = 0;
% F_right = (1/delta_x)*0.5*(F(i+1,j)+F(i,j));
% F_left =F_right;
% elseif i == N
% F_left = (1/delta_x)*0.5*(F(i,j)+F(i-1,j));
% % F_right = 0;
% F_right =F_left;
% else
% F_right = (1/delta_x)*0.5*(F(i+1,j)+F(i,j));
% F_left = (1/delta_x)*0.5*(F(i,j)+F(i-1,j));
% end
% 
% Unew(i,j) = U(i,j) - (delta_t)*(F_right-F_left)-(delta_t)*(G_up - ...
% G_down);
% end
% end
% U = Unew;
% surf(X,Y,Unew)
% pause(0.1);
% end

