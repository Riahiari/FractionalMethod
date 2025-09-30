function OneLayerParameters = FitOneLayerTPSFWavelengthVisual(Data, IRF, Guesses, times, PercentMode, startPercent, endPercent, lb, ub, rho, NumberOfFunctionCalls, functol, Plot)




NormalizeByArea = true;



DataSmoothed = smoothdata(Data, 2, 'movmean', 1);
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





    function Residual = ModelErrorOneLayer(parameters, Data, IRF, times, startIndex, endIndex, rho, NormalizeByArea)   

        %disp(parameters)

        mua1 = parameters(1); mus1 = parameters(2); A = parameters(3);

        Model = OneLayerTPSF(mua1, mus1, A, rho, times);
        Model(1) = 0;
        if any(IRF(:) ~= 0)
            Model = conv(Model, IRF, 'same');
        end
        Model = Model(startIndex:endIndex);

        if NormalizeByArea
            Model = Model ./ sum(Model);
            Data = Data ./ sum(Data);
        end

        Model = log(Model);
        Data = log(Data);

        lambda = 0;
        Factor = 1e0;
        
        Residual = norm(Model - Data, 'fro')^2 * Factor + lambda * norm(parameters, 'fro')^2;


        set(h2, 'YData', Model);
        hold on;
        set(h1, 'YData', Data);
        hold off;
        drawnow; 




    end








for i = 1:size(DataNormalized, 1)

    fprintf("Wavelength %i / %i\n", i, size(DataNormalized, 1))

    costFunction = @(params) ModelErrorOneLayer(params, fitData{i}, IRF(i,:), times, startPercentIndex(i), endPercentIndex(i), rho, NormalizeByArea);
    
    options = optimset('Display','None','MaxFunEvals',NumberOfFunctionCalls,'TolFun',functol,'TolX',functol,'MaxIter',NumberOfFunctionCalls,'UseParallel','Always'); % Options Parameter for homogenous fitting. Fitting is fast so it can have small TolX.

    [OptimizedParameters(i,:), FinalResidual] = fminsearchbnd(costFunction, Guesses(i,:), lb(i,:), ub(i,:), options);

end

OneLayerParameters = OptimizedParameters;





end