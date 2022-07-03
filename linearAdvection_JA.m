%=========================================================================
% Matlab program to solve the 1D Linear Advection equations using a choice
% of five explicit finite difference schemes. First Order Upwind,
% Lax-Friedrichs, Lax-Wendroff, Adams Average (Lax-Friedrichs) and Adams
% Average (Lax-Friedrichs). 
% Velocity is constant.
% A heuristic time step is used.
% Uses periodic boundary conditions (solutions reappears at the opposite
% side of the figure window).
% Initial conditions are discontinuous. 
% Plots the concentration (dependent variable) at each time step.
%
% Additional info
% The Adams Average scheme was devised by the author (James Adams) in 2014.
% The concept is trivial and has been applied to both the Lax-Friedrichs
% and Lax-Wendroff schemes. 
% Numerical experiments have shown that the Adams Average improves the
% performance of each scheme (given that the correct parameters are used).
% Recommended parameters (feel free to change these),
% Adams Average (Lax-Friedrichs) - A = 4 B = 4 C = 1
% Adams Average (Lax-Wendroff) - A = 1 B = 1 C = 5
% =========================================================================

function linearAdvection_JA
clear  % Clears workspace
clc   % Clears history
clf   % Clears figure window

N = 101;  % Defines the no. of grid points
p = 0;   % Sets the lower end of the domain
q = 100; % Sets the upper end of the domain
v = 0.8;  % defines velocity of the water (meters)
dx = (q - p)/(N - 1); % Calculates the grid spacing
dt = 0.8; % Defines the time step.
t = 0;  % Sets the initial time to zero
runtime = 100;  % Defines time program runs for
s = 0.4;    % Sets safety factor to 0.4 (Advised not to change this)
pausetime = 0.1;   % Pause time between animation

x = p : dx : q; % Calculates the value of x using the domain and the grid spacing

for i = 1 : N
    un(i) = U(x(i));  % Defines the initial pollutant using initial conditions
end

fprintf('Select a scheme from the list and enter the corrosponding number \n\nFirst Order Upwind - 1\nLax-Friedrichs - 2 \nLax-Wendroff - 3\nAdams Average (Lax-Friedrichs) - 4\nAdams Average (Lax-Wendroff) - 5\n\n')   
Choice = input('Scheme no. = ');   % Select a scheme to run

if Choice == 2 || Choice == 3   
    A = 1;
    B = 1;          % <<<< DO NOT CHANGE THESE VALUES 
    C = 0;
elseif Choice == 4 
    A = 2;
    B = 2;         % <<<<< Values can be changed
    C = 1;
elseif Choice == 5 
    A = 1;
    B = 1;        % <<<<<<<< Values can be changed 
    C = 5;
end
% ===================== Runs similation ===================================
while t < runtime  % Halts program when inequality is violated
    % Applies periodic boundary conditions depending on choice
    if Choice == 1
        boundary1 = un(N);
        boundary2 = un(1);
        un = [boundary1 un(1:N-1) boundary2];
    elseif Choice == 2 || Choice == 3 || Choice == 4 || Choice == 5
        boundary1 = un(N);
        boundary2 = un(1);
        un = [boundary1 un boundary2];
    end
    
    dt = s * dx / v;  % Calculates time step 
    t = t + dt;   % Adds time step to time 
    
    c = dt/dx;  % Calculates this here to provide efficiency
    
    if Choice == 1  % Applies First Order Upwind scheme
        for j = 2 : N + 1
            u(j) = un(j) - v * c * (un(j) - un(j-1)); 
        end
    elseif Choice == 2 || Choice == 4  % Applies either Lax-Friedrichs or Adams Average (Lax-Friedrichs) scheme
        for j = 2 : N + 1
            u(j) = (A*un(j+1) + B*un(j-1) + C*un(j))/(A+B+C) - v * c * (un(j+1) - un(j-1)); 
        end
    elseif Choice == 3 || Choice == 5  % Applies either Lax-Wendroff or Adams Average (Lax-Wendroff) scheme
        for j = 2 : N + 1
            u(j) = (A*un(j+1) + B*un(j-1) + C*un(j))/(A+B+C) - v * c * (un(j+1) - un(j-1))/2 + c^2 * v^2 * 0.5 * (un(j+1) - 2*un(j) + un(j-1)); 
        end
    end
    
    plot(x,u(2:N+1),'b')    % Plots pollutant (dependent variable)
    xlabel('x [m]')         % Adds appropiate labels
    ylabel('U [m]')
    
    % Adds title depending on scheme choice
    
    if Choice == 1
        title('Solves 1D Advection Equation using FOU Scheme','Fontsize',10)
    elseif Choice == 2
        title('Solves 1D Advection Equation using Lax-Friedrichs Scheme','Fontsize',10)
    elseif Choice == 3
        title('Solves 1D Advection Equation using Lax-Wendroff Scheme','Fontsize',10)
    elseif Choice == 4
        title('Solves 1D Advection Equation using Adams Average (Lax-Friedrichs) Scheme','Fontsize',10)
    elseif Choice == 5
        title('Solves 1D Advection Equation using Adams Average (Lax-Wendroff) Scheme','Fontsize',10)
    end
    
    axis([p q min(u) 1.1*max(u)])  % Sets axis 

    pause(pausetime)   % Shows results for each timestep
    
    un = u(2:N+1);   % Ommits boundary conditions ready for the next time step
  
end
end
% ======================= End of program ==================================

function initial = U(x) % Defines initial conditions.
pollutantstart1 = 20;   % Start of pollutant concentration 1
pollutantend1 = 40;     % End of pollutant concentration 1
pollutantstart2 = 60;   % Start of pollutant concentration 2
pollutantend2 = 80;     % End of pollutant concentration 2

pollutantfloor = 0.1;  % Defines bottom of pollutant
pollutantpeak = 0.8;   % Defines top of pollutant
initial = pollutantfloor;   % Sets up bottom of pollutant 

% ============ Determins the areas in which the pollutant peaks ==========

    if x >= pollutantstart1 && x <= pollutantend1
        initial = pollutantpeak;
    elseif x >= pollutantstart2 && x <= pollutantend2
        initial = pollutantpeak;
    end
end
