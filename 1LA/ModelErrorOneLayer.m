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
    
    end