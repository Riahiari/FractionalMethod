function TwoLayerParameters = FitTwoLayerTPSFWavelength(Data, Guesses, times, PercentMode, startPercent, endPercent, lb, ub, L, rho, NumberOfFunctionCalls, functol, PlotTrueFalse)

% Rho can be a vector



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

parfor i = 1:size(DataNormalized, 1)

    if PlotTrueFalse
    fprintf("Wavelength %i / %i\n", i, size(DataNormalized, 1))
    end

    costFunction = @(params) ModelErrorTwoLayer(params, fitData{i}, times, startPercentIndex(i), endPercentIndex(i), L, rho, NormalizeByArea);
    
    options = optimset('Display','None','MaxFunEvals',NumberOfFunctionCalls,'TolFun',functol,'TolX',0.1,'MaxIter',10000,'UseParallel','Always'); % Options Parameter for homogenous fitting. Fitting is fast so it can have small TolX.

    [OptimizedParameters(i,:), FinalResidual] = fminsearchbnd(costFunction, Guesses(i,:), lb(i,:), ub(i,:), options);

end

TwoLayerParameters = OptimizedParameters;



end