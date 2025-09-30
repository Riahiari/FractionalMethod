function TwoLayerParameters = FitTwoLayerTPSFWavelengthVisual(Data, Guesses, times, PercentMode, startPercent, endPercent, lb, ub, L, rho, NumberOfFunctionCalls, functol, Plot)






NormalizeByArea = true;



order = 3;
framelen = 7;

DataDenoised = zeros(size(Data));
DataSmoothed = zeros(size(Data));

for i=1:size(Data,1)
    DataDenoised(i,:) = Data(i,:);
    %DataDenoised(i,:) = denoise(Data(i,:));
    %DataDenoised( = wdenoise(Data(i,:));
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


waveIndex = 1;

figure(3);
h2 = plot(fitData{waveIndex}, 'r');
hold on;
h1 = plot(fitData{waveIndex}, 'b');
hold off;
drawnow;





    function Residual = ModelErrorTwoLayer(parameters, Data, times, startIndex, endIndex, L, rho, NormalizeByArea)   

        mua1 = parameters(1); mua2 = parameters(2); mus1 = parameters(3); mus2 = parameters(3); A = parameters(4);

        Model = TwoLayerTPSFifft(mua1, mua2, mus1, mus2, A, L, rho, times);
        Model = Model(startIndex:endIndex);

        if NormalizeByArea
            Model = Model ./ sum(Model);
            Data = Data ./ sum(Data);
        end

        Model = log(Model);
        Data =  log(Data);

        lambda = 0;
        Factor = 1e0;
        
        Residual = norm(Model - Data, 'fro')^2 * Factor + lambda * norm(parameters, 'fro')^2;


        set(h2, 'YData', real(Model));
        hold on;
        set(h1, 'YData', real(Data));
        hold off;
        drawnow; 



    end








for i = 1:size(DataNormalized, 1)

    fprintf("Wavelength %i / %i\n", i, size(DataNormalized, 1))

    costFunction = @(params) ModelErrorTwoLayer(params, fitData{i}, times, startPercentIndex(i), endPercentIndex(i), L, rho, NormalizeByArea);
    
    options = optimset('Display','None','MaxFunEvals',NumberOfFunctionCalls,'TolFun',functol,'TolX',functol,'MaxIter',NumberOfFunctionCalls,'UseParallel','Always'); % Options Parameter for homogenous fitting. Fitting is fast so it can have small TolX.

    [OptimizedParameters(i,:), FinalResidual] = fminsearchbnd(costFunction, Guesses(i,:), lb(i,:), ub(i,:), options);

end

TwoLayerParameters = OptimizedParameters;



end