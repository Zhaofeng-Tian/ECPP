function [fitresult, gof] = createFit(nnx, ylist1,SmoothingParam)
[xData, yData] = prepareCurveData( nnx, ylist1 );
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = SmoothingParam;
[fitresult, gof] = fit( xData, yData, ft, opts );


