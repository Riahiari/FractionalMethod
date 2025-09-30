function ScatteringCoeffs = DefineScatteringCoefficients(a, b, Wavelengths, lambda0)

ScatteringCoeffs = a.*(Wavelengths./lambda0).^(-b); % Mie scattering law

end