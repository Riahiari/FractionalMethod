function ScatterParameters = fitScatterModel(ScatteringCoeffs, wavelengths, lambda0, Plot)

if Plot
figure(5);
h1 = plot(wavelengths, log(ScatteringCoeffs), 'b');
hold on;
h2 = plot(wavelengths, log(ScatteringCoeffs), 'r');
end

% Linear Fitting

xValues = log(wavelengths./lambda0);
yValues = log(ScatteringCoeffs);

FitValues = polyfit(xValues, yValues, 1);

ScatterParameters = zeros(1, 2);

ScatterParameters(1) = exp(FitValues(2));
ScatterParameters(2) = -FitValues(1);

Model = DefineScatteringCoefficients(ScatterParameters(1), ScatterParameters(2), wavelengths, lambda0);

if Plot
set(h2, 'YData', log(Model));
end

hold off;

end