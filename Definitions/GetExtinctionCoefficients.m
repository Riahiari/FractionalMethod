function extinctions = GetExtinctionCoefficients(Wavelengths)  % Should later change to include variable chromophores.
   
    SPEdata = load('spe.mat'); % Load extinction coefficients.

    BoneAbsorptionData = load("mua_bone.mat");
    BoneMua = BoneAbsorptionData.BoneMua;

    spe = SPEdata.spe;
    
    wavelengthsRawSPE = spe(:,1);
    extinctionsRawSPE = [spe(:,4),spe(:,2),spe(:,3)]; % Extinction coefficients of interest [Water, Deoxy, Oxy]
    extinctionsSPE = interp1(wavelengthsRawSPE, extinctionsRawSPE, Wavelengths);

    wavelengthsRawBone = BoneMua(:,1);
    extinctionsRawBone = BoneMua(:,2);

    extinctionsBone = interp1(wavelengthsRawBone, extinctionsRawBone, Wavelengths)';

    extinctions = [extinctionsSPE]; % Combine extinction coefficients

end