
function AbsBPred = CalculateAbsB(AbsTPred, AbsRecPred, f, SmoothData)

%%%%%%%% SMOOTH DATA %%%%%%%%%%
if SmoothData
for u=1:size(AbsRecPred, 2)
    AbsTPred(:,u) = sgolayfilt(AbsTPred(:,u), 4, 31);
    AbsRecPred(:,u) = sgolayfilt(AbsRecPred(:,u), 4, 31);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% SOLVE FOR MU_A^B %%%
for j=1:size(AbsRecPred, 2)
AbsBPred(:,j) = (AbsRecPred(:,j) - f.*AbsTPred(:,j)) ./ (1 - f);
end
%%%%%%%%%%%%%%%%%%%%%%%%