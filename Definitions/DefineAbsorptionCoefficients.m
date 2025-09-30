function AbsorptionCoeffs = DefineAbsorptionCoefficients(c, extinctions)

    AbsorptionCoeffs = extinctions * c'; 

end