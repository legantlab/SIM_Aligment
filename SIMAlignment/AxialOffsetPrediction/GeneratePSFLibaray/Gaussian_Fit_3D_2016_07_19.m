function F = Gaussian_Fit_3D_2016_07_19(x,PeakData,dims)

F = (x(1) .* exp(-((PeakData(:,2) - x(2)).^2./(2*x(5)^2) + (PeakData(:,3) - x(3)).^2./(2*x(6)^2) + (PeakData(:,4) - x(4)).^2./(2*x(7)^2))));
F=reshape(F,dims(1),dims(2),dims(3));