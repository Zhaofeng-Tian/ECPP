function [x0,y0,x1,y1,x2,y2,x3,y3] = calculate(dx,x,y,theta)

w=0.4/dx;
L=0.8/dx;
lf=0/dx;
lr=0/dx;

kk=0.5;
p=x-L*kk*cos(theta);
q=y-L*kk* sin(theta);
x0=p-w/2* sin(theta)-lr*cos(theta);
y0=q+w/2*cos(theta)-lr*sin(theta);
x1=p+w/2* sin(theta)-lr*cos(theta);
y1=q-w/2 * cos(theta)-lr*sin(theta);

p=x+L*(1-kk)*cos(theta);
q=y+L*(1-kk)* sin(theta);
x2=p-w/2*sin(theta)+lf*cos(theta);
y2=q+w/2*cos(theta)+lf*sin(theta);
x3=p+w/2* sin(theta)+lf*cos(theta);
y3=q-w/2*cos(theta)+lf*sin(theta);

