
    function Residual = ModelErrorTwoLayer(parameters, Data, times, startIndex, endIndex, L, rho, NormalizeByArea)   

    % Rho can be a vector

        mua1 = parameters(1); mua2 = parameters(2); mus1 = parameters(3); mus2 = parameters(3); A = parameters(4);

        Model = TwoLayerTPSF(mua1, mua2, mus1, mus2, A, L, rho, times);
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

    end


