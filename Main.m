clear


addpath("Definitions/");
addpath("1LA/");
addpath("2LA/");
addpath("deconvolution/");
addpath("Fractional Data Loading/");


addpath("../LabComputerData/Fractional Slab Data Single Wavelength/Sep 3/");
addpath("../LabComputerData/Fractional Slab Data Single Wavelength/Sep 4/");
addpath("../LabComputerData/Fractional Slab Data Single Wavelength/Sep 10/");
addpath("../LabComputerData/Fractional Slab Data Single Wavelength/Sep 16/");
addpath("../LabComputerData/Fractional Slab Data Single Wavelength/Sep 18/");
addpath("../LabComputerData/Fractional Slab Data Single Wavelength/Sep 23/");
addpath("../LabComputerData/Fractional Slab Data Single Wavelength/Sep 29/");
addpath("../LabComputerData/Fractional Slab Data Single Wavelength/");


tstart = 0;
tend = 5e-9 * 1e9;
NumberOfTimeGates = 200;
time = linspace(tstart, tend, NumberOfTimeGates);

startPercent1 = 70;
endPercent1 = 1;
startPercent2 = -10;
endPercent2 =  1;

PercentMode = 1;

rhoValues = [5, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 23, 26, 30, 33, 36, 40];
rhoIndex1 = 5;
rhoIndex2 = 15;
sdd1 = rhoValues(rhoIndex1);
sdd2 = rhoValues(rhoIndex2);

UseMamadouDenoising = false;

    wavelengths = [701:2:851];
    %waveIndex = 1:76;
    %wavelengths = wavelengths(waveIndex);
    extinctions = GetExtinctionCoefficients(wavelengths);
    lambda0 = 800;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MULTI    WAVELENGTH     TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if true

    AriaDataStruct = load("RadialFluxDataOverAllWavelengths.mat");
    AriaData = AriaDataStruct.RadialFluxDataOverAllWavelengths;

    UseIRF = false;
    MamadouDenoising = false;

    if UseIRF
        IRF = zeros(size(squeeze(AriaData(1,:,:,1,1))))';
        IRF(:,1:15) = 1000;
    else
        IRF = squeeze(zeros(size(AriaData(1, :, :, 1, 1))))';
    end


%%% GUESSES AND LIMITS %%%

lbVals = [0.1, 10, 10,        0.5, 0.5, 1];
ubVals = [1.3, 90 ,90,        5, 5, 1];
guessVals = [0.1, 30, 30,     1, 1, 1];


GuessAndLimits = DefineGuessAndLimits(guessVals, lbVals, ubVals, extinctions, wavelengths, lambda0);
guess = GuessAndLimits(:,:,1);
lb = GuessAndLimits(:,:,2);
ub = GuessAndLimits(:,:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%% TOP LAYER GRADIENT %%%%%%%

DataStruct1 = load("RadialFluxDataOverAllWavelengths.mat");

    NumberOfAverages = 3;
    NumberOfConcs1 = 8;

for j=1:NumberOfAverages
for i=1:NumberOfConcs1
AbsTopIndex = i;
AbsBottomIndex = 4;

Data1 = LoadSourceDetectorDataForConcentrationIteration(DataStruct1, rhoIndex1, [AbsTopIndex, AbsBottomIndex, j], IRF, MamadouDenoising);
Data2 = LoadSourceDetectorDataForConcentrationIteration(DataStruct1, rhoIndex2, [AbsTopIndex, AbsBottomIndex, j], IRF, MamadouDenoising);


TopLayerFit = FitOneLayerTPSFWavelength(Data1, IRF, guess, time, PercentMode, startPercent1, endPercent1, lb, ub, sdd1, 1e10, 1e-10, 0);
BestGuessOneLayer = FitOneLayerTPSFWavelength(Data2, IRF, guess, time, PercentMode, startPercent2, endPercent2, lb, ub, sdd2, 1e10, 1e-10, 0);


AbsT(:,i,j) = TopLayerFit(:,1);
AbsRec(:,i,j) = BestGuessOneLayer(:,1);

%AbsT(:,i) = sgolayfilt(AbsT(:,i), 4, 25);
%AbsRec(:,i) = sgolayfilt(AbsRec(:,i), 4, 25);

end

AbsT = mean(AbsT, 3);
AbsRec = mean(AbsRec, 3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Calculate f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = CalculateFraction(AbsT, AbsRec, 1, true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(20)
clf
plot(AbsT(1,:),'r')
hold on
plot(AbsRec(1,:),'y')
hold off




%%%%%%%%%%% PREDICTION STAGE %%%%%%%%%%%

    RealWaterConcs = linspace(0.5,0.9,8);
    SkinSats = linspace(60,70,5);
    RealBrainSats = linspace(50,80,8);
    HbT = 55;

AbsRecPred = zeros(length(wavelengths), NumberOfConcs1, NumberOfAverages);
AbsTPred = zeros(size(AbsRecPred,1), NumberOfAverages);
RealTimes = 1:size(AbsRecPred,2);

DataStruct2 = load("RadialFluxDataOverAllWavelengths.mat");

for j=1:size(AbsRecPred,3)
for i=1:size(AbsRecPred,2)
AbsTopIndex = 1;
AbsBottomIndex = i;

Data1 = LoadSourceDetectorDataForConcentrationIteration(DataStruct2, rhoIndex1, [AbsTopIndex, AbsBottomIndex, j], IRF, MamadouDenoising);
Data2 = LoadSourceDetectorDataForConcentrationIteration(DataStruct2, rhoIndex2, [AbsTopIndex, AbsBottomIndex, j], IRF, MamadouDenoising);


TopLayerFit = FitOneLayerTPSFWavelength(Data1, IRF, guess, time, PercentMode, startPercent1, endPercent1, lb, ub, sdd1, 1e10, 1e-10, 0);
BestGuessOneLayer = FitOneLayerTPSFWavelength(Data2, IRF, guess, time, PercentMode, startPercent2, endPercent2, lb, ub, sdd2, 1e10, 1e-10, 0);


AbsTPred(:,i,j) = TopLayerFit(:,1);
AbsRecPred(:,i,j) = BestGuessOneLayer(:,1);
ScatRecPred(:,i,j) = BestGuessOneLayer(:,2);

AbsBPredUnsmooth(:,i,j) = CalculateAbsB(AbsTPred(:,i,j), AbsRecPred(:,i,j), f, false);
AbsBPred(:,i,j) = CalculateAbsB(AbsTPred(:,i,j), AbsRecPred(:,i,j), f, true);
end


end

AbsBPred = mean(AbsBPred, 3);
AbsBPredUnsmooth = mean(AbsBPredUnsmooth, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%% B Absorption ANALYSIS %%%%%

for i=1:size(AbsRecPred, 2)
AbsRealB(:,i) = DefineAbsorptionCoefficients([RealWaterConcs(1),HbT*(1-RealBrainSats(i)/100),HbT*RealBrainSats(i)/100], extinctions);
end

ErrorArray = abs(abs(AbsRealB) - abs(AbsBPred)) ./ AbsRealB .* 100;


figure(3)
clf
plot(RealTimes, AbsBPred(1,:), 'r')
hold on
plot(RealTimes, AbsRealB(1,:), 'b')
hold off

figure(4)
clf
plot(AbsBPred(:,1), 'r')
hold on
plot(AbsRealB(:,1), 'b')
hold on
plot(AbsBPredUnsmooth(:,1), 'y')
hold off

disp("MEAN ERRORS AT EACH SIMULATION ITERATION")
disp(mean(ErrorArray))


end





%%%%%%%%%%%%% CONCENTRATION DETERMINATIONS %%%%%%%%%%%%%%%

for i = 1:size(AbsBPred,2)

    ConcValues(i,:) = GetConcentrationsFromAbsorbance(AbsBPred(:,i), extinctions, [0.5, 15, 15], [0.2, 3, 3], [3.5, 60, 60], 1, wavelengths, 1, 0);

end

disp(ConcValues)


PredBrainWater = ConcValues(:,1);
PredBrainSats = ConcValues(:,3) ./ (ConcValues(:,2) + ConcValues(:,3)) .* 100;
BrainWaterError = mean(abs(PredBrainWater' - RealWaterConcs) ./ RealWaterConcs .* 100);
BrainSatsError = mean(abs(PredBrainSats' - RealBrainSats) ./ RealBrainSats .* 100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






figure(7)
clf
plot(PredBrainWater, 'r')
hold on
plot(RealWaterConcs, 'b')
hold off


figure(8)
clf
plot(PredBrainSats, 'r')
hold on
plot(RealBrainSats, 'b')
hold off



if true
    figure(12)
    clf
    plot(AbsBPred(:,1),'r')
    hold on
    plot(AbsRealB(:,1), 'b')
    hold on
    plot(AbsBPredUnsmooth(:,1), 'y')
    hold off
end




if true

for j=1:NumberOfConcs1
for q=1:NumberOfAverages
AbsTopIndex = 1;
AbsBottomIndex = j;
LValues = [10,11,12,13,14,15,12,14];
Data1 = LoadSourceDetectorDataForConcentrationIteration(DataStruct2, rhoIndex1, [AbsTopIndex, AbsBottomIndex, q], IRF, MamadouDenoising);
Data2 = LoadSourceDetectorDataForConcentrationIteration(DataStruct2, rhoIndex2, [AbsTopIndex, AbsBottomIndex, q], IRF, MamadouDenoising);
TwoLayerModelAbsorptions(:,j,q) = BottomLayerAbsorptionFromTwoLayerModel(Data2, AbsTPred(:,j,q), ScatRecPred(:,j,q), time, -10, 1, LValues(j), sdd2, extinctions, wavelengths, lambda0);
end
end

ErrorArrayTwoLayerModel = abs(abs(AbsRealB) - abs(mean(TwoLayerModelAbsorptions,3))) ./ AbsRealB .* 100;
disp("Two Layer Model Errors")
disp(mean(ErrorArrayTwoLayerModel))

TLMAM = mean(TwoLayerModelAbsorptions, 3);
Error3 = mean(abs(abs(TLMAM) - abs(AbsRealB(:,1))) ./ AbsRealB(:,1) .* 100);
disp("Averaged Error")
disp(Error3)


if true
    figure(23)
    clf
    plot(AbsBPred(:,1),'r')
    hold on
    plot(AbsRealB(:,1), 'b')
    hold on
    plot(AbsBPredUnsmooth(:,1), 'y')
    hold on
    plot(TLMAM, 'g')
    hold off
end


end

