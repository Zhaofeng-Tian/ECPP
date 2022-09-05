%% 
% Author: Zhaofeng Tian
% Paper: Cutting Edge Method to Cut the Edge for Robot Mowing
% Email: Shoguntian@gmail.com

%%
close all
clear all
%% ***********************Parameter Initializatin

L=0.8;% length
B=0.4;% width
Rmax=sqrt(L^2+B^2)/2;% R1 big disk
RmGE=0.15;% R3 mowing disk
RmB=B/2; % R2 small disk

dx=0.01;% map resolution x
dy=0.01;% map resolution y

load('datarenyi.mat');
load('datarenyi2.mat');

d6=[Sxx' Syy'];
dddd=1000;
ph=0.99;
[fitresult, gof] =smoothing(d6,dddd,ph) ;
Sx=0:dx:10;
Sy=fitresult(Sx);

figure(6),
plot(Sxx,Syy);hold on
plot(Sx,Sy);

Nx=ceil(max(Sx)/dx);
Ny=ceil(max(abs(Sy))/dy)*4;
starty=ceil(max(abs(Sy))/dy)*2;

A=zeros(Ny,Nx);

for i=1:Nx
    nnx(i)=i;
    nny(i)=starty+ceil(Sy(i)/dy);
    A(nny(i),i)=1;
    
    if i>1
        a1=min(nny(i-1),nny(i));
        a2=max(nny(i-1),nny(i));
        
        A(a1:a2,i-1)=1;
    end
end
%***********SE
R=ceil(Rmax/dx);
SE1=strel('disk',R,0)  
figure(3),imagesc(SE1.Neighborhood);axis equal 
% %% ****************Closing Operation
A_close=imclose(A,SE1); % dilate
%% *********************** Boundary Preprocessing
figure(1),
x=nnx;
y=nny;
plot(x,y,'Color',[0.00,0.78,0.55],'linewidth',5);

hold on;axis equal 
for i=1:Nx
    [a ,b]=find(A_close(:,i)==1);
    ylist_low1(i)=a(1);
    ylist_up1(i)=a(end);
end
%y=ylist_up1;
%plot(nnx,ylist_up1,'+','Color',[0.89,0.81,0.34]);


%hold on;
y=ylist_low1;
plot(nnx,ylist_low1,'o','Color',[0.63,0.32,0.18]);
xlim([400 600])
ylim([400 600])
hold on;axis equal 
%% ****************************************** Imaginary Mowing Disk Planning
R=ceil(RmGE/dx);
SE0=strel('disk',R,0)

A_close_dilate=imdilate(A_close,SE0); 
figure(1),
for i=1:Nx
    [a ,b]=find(A_close_dilate(:,i)==1);
    ylist_low0(i)=a(1);
end
plot(nnx,ylist_low0,'k--');hold on;axis equal

%% ********************************************** "Dilation-by-circumcircle" - "Big disk" planning
A_close_dilate=imdilate(A_close,SE1); 
figure(2),

subplot(131),imagesc(A);axis equal
set(gca,'YDir','normal');
ylim([200,800]) 
 
 
subplot(132),imagesc(A_close);axis equal 
set(gca,'YDir','normal');
ylim([200,800]) 
subplot(133),imagesc(A_close_dilate);axis equal 
set(gca,'YDir','normal');
ylim([200,800]) 


figure(1),
for i=1:Nx
    [a ,b]=find(A_close_dilate(:,i)==1);
    ylist_low2(i)=a(1);
end
plot(nnx,ylist_low2,'b--');hold on;axis equal 
%% ****************************************** "Small disk" Planning
R=ceil(RmB/dx);
SE2=strel('disk',R,0)  
A_close_dilate=imdilate(A_close,SE2); 
figure(1),
for i=1:Nx
    [a ,b]=find(A_close_dilate(:,i)==1);
    ylist_low3(i)=a(1);
end
plot(nnx,ylist_low3,'g--');hold on;axis equal 
%% **************************************** Convexity Calculation

figure,
subplot(311)
%plot(nnx,ylist_low1,'r');hold on;
[fitresult, gof] = createFit(nnx, ylist_low1,0.001);
ysmooth=fitresult(nnx);ysmooth=ysmooth(:);

plot(nnx,ysmooth,'b');axis equal;hold on;
%*********** First order
tempy=[ysmooth(2)-ysmooth(1);ysmooth(2:end)-ysmooth(1:end-1)];
dxx=1;
dy1=tempy/dxx;

subplot(312)
plot(nnx,dy1,'b');hold on;
plot(nnx,dy1*0,'k--');hold on;legend('First-order derivative')
%*********** Second order

tempy=[dy1(2)-dy1(1);dy1(2:end)-dy1(1:end-1)];
dxx=1;
dy2=tempy/dxx;
subplot(313)
plot(nnx,dy2,'b');hold on;
plot(nnx,dy2*0,'k--');hold on;legend('Second-order derivative')
%% ************************************* "Big and small disk" Planning
for i=1:Nx
    if dy2(i)>0
        ylist_low_max_min(i)=ylist_low3(i);  % small disk
    else
        ylist_low_max_min(i)=ylist_low2(i);  % big disk
    end
end
figure(1)
plot(nnx,ylist_low_max_min,'m-');hold on;axis equal %

[fitresult, gof] = createFit(nnx, ylist_low_max_min,0.2);
ylist_low_max_min_smooth=fitresult(nnx);ylist_low_max_min_smooth=ylist_low_max_min_smooth(:);

plot(nnx,ylist_low_max_min_smooth,'c-');hold on;axis equal %

%% *************************************** " Sliding Chopstick " Planning
Length_L=L/dx;  
Width=RmB/dx;
for i=1:Nx-Length_L
    for j=i+1:i+Length_L
        Lnew(j-i)=sqrt((ylist_low1(j)-ylist_low1(i))^2+(nnx(i)-nnx(j))^2);
    end
    temp=abs( Lnew-Length_L);% Find closest point that makes the distance equal to L
    [aa ,bb]=min(temp);
    x0=(i+i+bb)/2;
    y0=(ylist_low1(i)+ylist_low1(i+bb))/2;
    kab=(ylist_low1(i+bb)-ylist_low1(i))/bb;
    kabchui=-1/kab;
    
    if kab>0
        xnew(i)=x0+Width/sqrt(kabchui^2+1);
        ynew(i)=y0-Width/sqrt(kabchui^2+1)*abs(kabchui);
    elseif kab==0
        xnew(i)=x0;
        ynew(i)=y0-Width;
    else
        xnew(i)=x0-Width/sqrt(kabchui^2+1);
        ynew(i)=y0-Width/sqrt(kabchui^2+1)*abs(kabchui);
    end
    
    if 1
        x0=i+1;
        y0=(ylist_low1(i)+ylist_low1(i+2))/2;
        %         kab=(ylist_low1(i+2)-ylist_low1(i))/2;
        kab=dy1(i);
        kabchui=-1/kab;
        if kab>0
            xnew2(i)=x0+Width/sqrt(kabchui^2+1);
            ynew2(i)=y0-Width/sqrt(kabchui^2+1)*abs(kabchui);
        elseif kab==0
            xnew2(i)=x0;
            ynew2(i)=y0-Width;
        else
            xnew2(i)=x0-Width/sqrt(kabchui^2+1);
            ynew2(i)=y0-Width/sqrt(kabchui^2+1)*abs(kabchui);
        end
    end
end
%% 



Nmin=max([xnew(1) xnew2(1)]);
Nmax=min([xnew(end) xnew2(end)]);
Nlist=ceil(Nmin):1:fix(Nmax);
Y=[interp1(xnew,ynew,Nlist)
    interp1(xnew2,ynew2,Nlist)];
Yend=min(Y);
Yend=[ylist_low_max_min_smooth(1:Nlist(1)-1)' Yend ylist_low_max_min_smooth(Nlist(end)+1:max(nnx))'];
Nlist=nnx;
[fitresult, gof] = createFit(Nlist, Yend,0.2);
ylist_line_smooth=fitresult(Nlist);

% figure
% plot(xnew,ynew,'k-');hold on;axis equal %
% plot(xnew2,ynew2,'m-');hold on;axis equal %
% plot(Nlist,Yend,'g-');hold on;axis equal %
% plot(Nlist,ylist_line_smooth,'c-');hold on;axis equal %

figure(1)
plot(Nlist,Yend,'r-');hold on;axis equal %

%% *************************** Data Recording
figure
plot(nnx,nny,'k-');hold on;axis equal % raw boundary
xlim([0,1000])
ylim([300,700])
plot(nnx,ylist_low0,'k--');hold on;axis equal % imaginary mowing disk planning
xlim([0,1000])
ylim([300,700])
plot(nnx,ylist_low1,'b--');hold on;axis equal % preprocessed boundary
xlim([0,1000])
ylim([300,700])

plot(nnx,ylist_low2,'r--');hold on;axis equal % big disk planning
xlim([0,1000])
ylim([300,700])
plot(nnx,ylist_low3,'g--');hold on;axis equal % small disk planning
xlim([0,1000])
ylim([300,700])
plot(nnx,ylist_low_max_min_smooth,'c-');hold on;% big and small disk planning
xlim([0,1000])
ylim([300,700])
% plot(nnx,ylist_low_max_min,'k-');hold on;axis equal 
plot(nnx,ylist_line_smooth,'m-');hold on;axis equal % sliding chopstick planning
xlim([0,1000])
ylim([300,700])

legend('Original boundary','Imaginary mowing disk','Preprocessed boundary','Big disk planning','Small disk planning','Big and small disk planning','Sliding chopstick planning')
%% ****************************** Uncut Area Calculation
figure
canyu.Rmax=0;
canyu.Rmin=0;
canyu.RmaxRmin=0;
canyu.line=0;

for i=1:length(nnx)
    canyu.Rmin=(ylist_low3(i)-ylist_low0(i))*1+canyu.Rmin;
    canyu.Rmax=(ylist_low2(i)-ylist_low0(i))*1+canyu.Rmax;
    canyu.RmaxRmin=(ylist_low_max_min_smooth(i)-ylist_low0(i))*1+canyu.RmaxRmin;
    canyu.line=(ylist_line_smooth(i)-ylist_low0(i))*1+ canyu.line;
end
data=-[canyu.Rmax canyu.Rmin canyu.RmaxRmin canyu.line];
data=data*(dx*dy);
c = categorical({'Big disk','Small disk','Big and small disk','Sliding chopstick'});
bar(data);

set(gca,'xticklabel',{'Big disk','Small disk','Big and small disk','Sliding chopsitck'});
xlabel('Planning Method')
ylabel('Uncut area [m/s^2]')


%% ***************************** Plot Trajectories
dx=dx*1.3;

figure
subplot(311)

[fitresult, gof] = createFit(nnx, ylist_low2,0.05);
 plot(nnx,nny,'k-');hold on;


plot(nnx,ylist_low1,'b--');hold on;axis equal


plot(nnx,ylist_low2,'r--');   hold on;axis equal 
plot(nnx,fitresult(nnx),'b-');hold on;axis equal


for i=1:1:max(nnx)-1
    KKK=fitresult(i+5)-fitresult(i);
    theta(i)=atan(KKK);
    x=i;
    y=fitresult(i);
    [x0(i),y0(i),x1(i),y1(i),x2(i),y2(i),x3(i),y3(i)] = calculate(dx ,x,y,theta(i));
    plot([x0(i) x1(i) x3(i) x2(i) x0(i)],[y0(i) y1(i) y3(i) y2(i) y0(i)],'b-');hold on;

end


subplot(312)

[fitresult, gof] = createFit(nnx, ylist_low_max_min_smooth,0.05);
plot(nnx,nny,'k-');hold on;

plot(nnx,ylist_low1,'b--');hold on;axis equal 


plot(nnx,ylist_low_max_min_smooth,'r--');   hold on;axis equal 


plot(nnx,fitresult(nnx),'b-');hold on;axis equal 


for i=1:1:max(nnx)-1
    KKK=fitresult(i+5)-fitresult(i);
    theta(i)=atan(KKK);
    x=i;
    y=fitresult(i);
    [x0(i),y0(i),x1(i),y1(i),x2(i),y2(i),x3(i),y3(i)] = calculate(dx ,x,y,theta(i));
    plot([i i+1 ],[fitresult(i),fitresult(i+1)],'k*')
    plot([x0(i) x1(i) x3(i) x2(i) x0(i)],[y0(i) y1(i) y3(i) y2(i) y0(i)],'b-');hold on;
end


subplot(313)

[fitresult, gof] = createFit(nnx, ylist_line_smooth,0.05);
plot(nnx,nny,'k-');hold on;


plot(nnx,ylist_low1,'b--');hold on;axis equal 


plot(nnx,ylist_line_smooth,'r--');   hold on;axis equal 


plot(nnx,fitresult(nnx),'b-');hold on;axis equal 

for i=1:3:max(nnx)-1
    KKK=fitresult(i+6)-fitresult(i);
    theta(i)=atan(KKK);
    x=i;
    y=fitresult(i);
    [x0(i),y0(i),x1(i),y1(i),x2(i),y2(i),x3(i),y3(i)] = calculate(dx ,x,y,theta(i));
    plot([x0(i) x1(i) x3(i) x2(i) x0(i)],[y0(i) y1(i) y3(i) y2(i) y0(i)],'b-');hold on;
end
%title

