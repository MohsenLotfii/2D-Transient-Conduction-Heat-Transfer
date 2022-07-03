
function heatNeumann(t_i,t_f,a,b,dx,dt,alpha)

% dt: step size in time dimension
% dx: step size in x axes
% a: left point of domain x
% b: right point of domain x
% alpha: take 1
% t_i: t initial, t_f: t final

n=(t_f-t_i)/dt;
m=(b-a)/dx;
lambda=alpha*dt/dx^2;
T=zeros(m+1,n+1);

x=a:dx:b;
t=t_i:dt:t_f;

u0 = x.^2+1+cos(pi.*x);
T(:,1)=u0; %initial value
u_left = T(2,1); %boundary cond.
u_right = T(m,1)+4*dx; %boundary cond.

for j=1:n
    T(1,j+1)=T(1,j)+lambda*(T(2,j)-2*T(1,j)+u_left); %boundary cond.
    for i=2:m
        T(i,j+1)=T(i,j)+lambda*(T(i+1,j)-2*T(i,j)+T(i-1,j));
    end
    T(m+1,j+1)=T(m+1,j)+lambda*(u_right-2*T(m+1,j)+T(m,j)); %boundary cond.
    u_left = T(2,j+1);
    u_right = T(m,j+1);
end

%exact soln

[xx,tt]=meshgrid(x,t);
exact=2.*tt+xx.^2+1+exp(-pi^2.*tt).*cos(pi.*xx);

[T(:,end) exact(end,:)']
