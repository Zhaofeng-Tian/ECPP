function [fitresult, gof] =smoothing(d6,dddd,ph) 
xp=d6(:,1)';yp=d6(:,2)';
[tmp I] = unique(d6(:,1), 'first');d6 = d6(I,:);xp=d6(:,1);yp=d6(:,2);
xxx=linspace(min(xp),max(xp),dddd);
yyy=interp1(xp,yp,xxx,'linear');
[xData, yData] = prepareCurveData( xxx, yyy );
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = ph;  
[fitresult, gof] = fit( xData, yData, ft, opts );

