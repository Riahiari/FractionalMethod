function BottomAbsorption = BottomLayerAbsorptionFromTwoLayerModel(Data, AbsT, ScatteringCoeff, time, startPercent, endPercent, SkullThickness, SourceDetectorDistance, extinctions, wavelengths, lambda0)


%% TOP LAYER FITTING %%


%% -------- SECOND LAYER INITIALIZATION WITH OPTIMIZATION ----------- %%

% First Round

if true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberOfChromophores = 3;

GuessVals =      [1,15,15,                 0.8,30,30,                    1,1, 1];
UpperBoundVals = [1,15,15,                 1.6, 150, 150,                  1,1, 1];
LowerBoundVals = [1,15,15,                 0, 0, 0,                      1,1, 1];

GuessAndLimits = DefineGuessAndLimits(GuessVals, LowerBoundVals, UpperBoundVals, extinctions, wavelengths, lambda0);

Guess = GuessAndLimits(:,:,1);
Guess(:,1) = AbsT;
Guess(:,3) = ScatteringCoeff;

LowerBound = GuessAndLimits(:,:,2);
LowerBound(:,1) = AbsT;
LowerBound(:,3) = ScatteringCoeff;

UpperBound = GuessAndLimits(:,:,3);
UpperBound(:,1) = AbsT;
UpperBound(:,3) = ScatteringCoeff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavetic = tic;
RelevantParameters = FitTwoLayerTPSFWavelength(Data, Guess, time, 1, startPercent, endPercent, LowerBound, UpperBound, SkullThickness, SourceDetectorDistance, lambda0, 0.000001, 0);
wavetoc = toc(wavetic);
 
fprintf("Fitting Time: %.2f\n", wavetoc)

AbsorptionsSecondLayer = RelevantParameters(:, 2);

Concs2 = GetConcentrationsFromAbsorbance(AbsorptionsSecondLayer, extinctions, GuessVals(NumberOfChromophores+1:2*NumberOfChromophores), LowerBoundVals(NumberOfChromophores+1:2*NumberOfChromophores), UpperBoundVals(NumberOfChromophores+1:2*NumberOfChromophores), 1, wavelengths, 0, 0);

BestGuessParameters = Concs2;

disp("BEST GUESS FROM TWO LAYERS")
disp(BestGuessParameters)
disp("OXYGENATION")
disp( Concs2(3) * 100 / ( Concs2(2) + Concs2(3) ) )

BottomAbsorption = AbsorptionsSecondLayer;
end

