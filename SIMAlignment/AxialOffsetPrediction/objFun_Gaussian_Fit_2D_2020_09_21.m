function F = objFun_Gaussian_Fit_2D_2020_09_21(x,PeakData)

F = sum(((PeakData(:,1) - (x(1) .* exp(-((PeakData(:,2) - x(2)).^2./(2*x(4)^2) + (PeakData(:,3) - x(3)).^2./(2*x(5)^2))))).^2));