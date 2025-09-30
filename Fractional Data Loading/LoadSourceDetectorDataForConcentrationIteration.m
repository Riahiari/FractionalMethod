function Data = LoadSourceDetectorDataForConcentrationIteration(DataStruct, RhoIndex, Indices, IRF, Denoising)

    % Load in Data %
    AriaDataStruct = DataStruct;
    AriaData = AriaDataStruct.RadialFluxDataOverAllWavelengths;

    % Apply Correct Concentration Indices based on Loaded Data Dimensions %
    idx = repmat({':'}, 1, ndims(AriaData));
    idx{1} = RhoIndex;
    for k = 1:numel(Indices)
        idx{k+3} = Indices(k);
    end
    
    Data = squeeze(AriaData(idx{:}))';

    % Convolve IRF %
    if any(IRF) ~= 0
        Data = conv(Data, IRF, 'same');
    end

    % Denoise Data %
    if Denoising
        for u=1:size(Data,1)
        Data(u,:) = denoised(Data(u,:));
        end
    end

end