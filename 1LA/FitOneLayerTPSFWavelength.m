function OneLayerParameters = FitOneLayerTPSFWavelength(Data, IRF, Guesses, times, PercentMode, startPercent, endPercent, lb, ub, rho, NumberOfFunctionCalls, functol, Plot)


% lbGlobal and ubGlobal are the upper and lower bounds on 'global'
% variables, such as a,b in mus = a(lambda/lambda0)^(-b), and the
% concentrations that determine mua. 


NormalizeByArea = true;

order = 3;
framelen = 7;

DataDenoised = zeros(size(Data));
DataSmoothed = zeros(size(Data));

for i=1:size(Data,1)
    DataDenoised(i,:) = Data(i,:);
    %DataDenoised(i,:) = denoise_only(Data(i,:));
    %DataDenoised(i,:) = wdenoise(Data(i,:));
    DataSmoothed(i,:) = DataDenoised(i,:);
    %DataSmoothed(i,:) = sgolayfilt(DataDenoised(i,:), order, framelen);
end
NormalizationFactor = max(DataSmoothed, [], 2);
DataNormalized = DataSmoothed ./ NormalizationFactor;

fitData = cell(size(DataNormalized, 1), 1); % Pre-allocate cell array

startPercentIndex = zeros(size(Data, 2), 1);
endPercentIndex = zeros(size(Data, 2), 1);


for i = 1:size(DataNormalized, 1)
    if PercentMode
    if startPercent <= 0
    startPercentIndex(i) = find( DataNormalized(i,:) > (-startPercent/100), 1, 'first' );
    elseif startPercent > 0
    startPercentIndex(i) = find( DataNormalized(i,:) > (startPercent/100), 1, 'last' );
    end
    endPercentIndex(i) = find( DataNormalized(i,:) > (endPercent/100), 1, 'last' );
    else
    startPercentIndex(i) = startPercent;
    endPercentIndex(i) = endPercent;
    end

    fitData{i} = DataSmoothed(i, startPercentIndex(i):endPercentIndex(i));
end

if Plot
waveIndex = 1;

figure(5);
h2 = plot(log(fitData{waveIndex}), 'r');
hold on;
h1 = plot(log(fitData{waveIndex}), 'b');
hold off;
drawnow;
end



parfor i = 1:size(DataNormalized, 1)

    if Plot
    fprintf("Wavelength %i / %i\n", i, size(DataNormalized, 1))
    end

    costFunction = @(params) ModelErrorOneLayer(params, fitData{i}, IRF(i,:), times, startPercentIndex(i), endPercentIndex(i), rho, NormalizeByArea);
    
    options = optimset('Display','None','MaxFunEvals',NumberOfFunctionCalls,'TolFun',functol,'TolX',functol,'MaxIter',NumberOfFunctionCalls,'UseParallel','Always'); % Options Parameter for homogenous fitting. Fitting is fast so it can have small TolX.

    [OptimizedParameters(i,:), FinalResidual] = fminsearchbnd(costFunction, Guesses(i,:), lb(i,:), ub(i,:), options);

end

OneLayerParameters = OptimizedParameters;





end