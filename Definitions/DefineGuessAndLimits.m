function GuessAndLimitsArray = DefineGuessAndLimits(GuessVals, LowerBoundVals, UpperBoundVals, extinctions, wavelengths, lambda0)

if length(GuessVals) == 6 || length(GuessVals) == 7

    NumberOfChromophores = length(GuessVals) - 3;

    GuessAbsorption1 = DefineAbsorptionCoefficients(GuessVals(1:NumberOfChromophores), extinctions);
    GuessScatter1 = DefineScatteringCoefficients(GuessVals(end-2), GuessVals(end-1), wavelengths, lambda0)';
    GuessAmplitude = GuessVals(end) * ones(1, length(wavelengths))';

    Guess = [GuessAbsorption1, GuessScatter1, GuessAmplitude];

    LowerAbsorption1 = DefineAbsorptionCoefficients(LowerBoundVals(1:NumberOfChromophores), extinctions);
    LowerScatter1 = DefineScatteringCoefficients(LowerBoundVals(end-2), LowerBoundVals(end-1), wavelengths, lambda0)';
    LowerAmplitude = LowerBoundVals(end) * ones(1, length(wavelengths))';
    
    LowerBound = [LowerAbsorption1, LowerScatter1, LowerAmplitude];
    
    UpperAbsorption1 = DefineAbsorptionCoefficients(UpperBoundVals(1:NumberOfChromophores), extinctions);
    UpperScatter1 = DefineScatteringCoefficients(UpperBoundVals(end-2), UpperBoundVals(end-1), wavelengths, lambda0)';
    UpperAmplitude = UpperBoundVals(end) * ones(1, length(wavelengths))';
    
    UpperBound = [UpperAbsorption1, UpperScatter1, UpperAmplitude];
    
    GuessAndLimitsArray = zeros(length(wavelengths), 3, 3);
    
    GuessAndLimitsArray(:,:,1) = Guess;
    GuessAndLimitsArray(:,:,2) = LowerBound;
    GuessAndLimitsArray(:,:,3) = UpperBound;


elseif length(GuessVals) == 9 || length(GuessVals) == 11

    NumberOfChromophores = (length(GuessVals) - 3)/2;

    GuessAbsorption1 = DefineAbsorptionCoefficients(GuessVals(1:NumberOfChromophores), extinctions);
    GuessAbsorption2 = DefineAbsorptionCoefficients(GuessVals(NumberOfChromophores+1:2*NumberOfChromophores), extinctions);
    GuessScatter1 = DefineScatteringCoefficients(GuessVals(end-2), GuessVals(end-1), wavelengths, lambda0)';
    GuessScatter2 = GuessScatter1;
    GuessAmplitude = GuessVals(end) * ones(1, length(wavelengths))';
    
    
    Guess = [GuessAbsorption1, GuessAbsorption2, GuessScatter1, GuessAmplitude];
    
    LowerAbsorption1 = DefineAbsorptionCoefficients(LowerBoundVals(1:NumberOfChromophores), extinctions);
    LowerAbsorption2 = DefineAbsorptionCoefficients(LowerBoundVals(NumberOfChromophores+1:2*NumberOfChromophores), extinctions);
    LowerScatter1 = DefineScatteringCoefficients(LowerBoundVals(end-2), LowerBoundVals(end-1), wavelengths, lambda0)';
    LowerScatter2 = LowerScatter1;
    LowerAmplitude = LowerBoundVals(end) * ones(1, length(wavelengths))';
    
    LowerBound = [LowerAbsorption1, LowerAbsorption2, LowerScatter1, LowerAmplitude];
    
    UpperAbsorption1 = DefineAbsorptionCoefficients(UpperBoundVals(1:NumberOfChromophores), extinctions);
    UpperAbsorption2 = DefineAbsorptionCoefficients(UpperBoundVals(NumberOfChromophores+1:2*NumberOfChromophores), extinctions);
    UpperScatter1 = DefineScatteringCoefficients(UpperBoundVals(end-2), UpperBoundVals(end-1), wavelengths, lambda0)';
    UpperScatter2 = UpperScatter1;
    UpperAmplitude = UpperBoundVals(end) * ones(1, length(wavelengths))';
    
    UpperBound = [UpperAbsorption1, UpperAbsorption2, UpperScatter1, UpperAmplitude];
    
    GuessAndLimitsArray = zeros(length(wavelengths), 4, 3);
    
    GuessAndLimitsArray(:,:,1) = Guess;
    GuessAndLimitsArray(:,:,2) = LowerBound;
    GuessAndLimitsArray(:,:,3) = UpperBound;

else
    
    disp(length(GuessVals))
    error('Invalid number of GuessVals. Expected 7 or 11 elements.');

end

end

