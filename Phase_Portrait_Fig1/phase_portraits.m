clear all
clc
% clf
n    = 2;
beta = 0.15;
K    = 10;
gamma= 10;
t = 0;
T = 150;
dt = 22; 
f = @(t,Y) [(beta*K^n)/(K^n+Y(2)^n)-(1/(60*gamma))*Y(1); (beta*K^n)/(K^n+Y(1)^n)-(1/(60*gamma))*Y(2)]; %Differential Equations to Update TFs
y1 = linspace(t,T,dt);
y2 = linspace(t,T,dt);
[x,y] = meshgrid(y1,y2);
size(x)
size(y)
u = zeros(size(x));
v = zeros(size(x));
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end
hh = quiver(x,y,u,v,2,'b'); figure(gcf)
ax = gca; ax.LineWidth = 1.1;
xlabel('Phenotype A','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Phenotype B','FontSize',16,'FontWeight','bold','Color','k')
xticks([t:75:T])
yticks([t:75:T])
axis tight equal;
circle(t,75,15,1.5,2.5,T) 
circle(75,t,15,0  ,1  ,T) 
hold on
plot(1:T,'k','LineWidth',2) 
hold off
function h = circle(x,y,r,sp,ep,nsegments)  
if nargin<4
    nsegments=50;
end
hold on
th = sp*pi:2*pi/nsegments:(ep*pi);
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,'r','LineWidth',2);
fill(xunit,yunit,'r')
colormap([0 0 0; 1 1 1]);
hold off
end