function Concentrations = GetConcentrationsFromAbsorbance(AbsorbanceVector, extinctionsUnsmoothed, ConcGuess, LowerBound, UpperBound, RegularFitting, wavelengths, DerivativeOrder, TrueValue)


if ~RegularFitting

order = 3;
frameLen = 15;

for i = 1:size(extinctionsUnsmoothed, 2)
    extinctions(:,i) = sgolayfilt(extinctionsUnsmoothed(:,i), order, frameLen);                      
    D1extinctionsUnsmoothed(:,i) = gradient(sgolayfilt(extinctionsUnsmoothed(:,i), order, frameLen));  
    D1extinctions(:,i) = sgolayfilt(D1extinctionsUnsmoothed(:,i), order, frameLen); 
    D2extinctions(:,i) = sgolayfilt(gradient(D1extinctions(:,i)), order, frameLen);
end

MuaDenoised = denoise_only(AbsorbanceVector);
%MuaDenoised = AbsorbanceVector;
Mua = sgolayfilt(MuaDenoised, order, frameLen);
D1MuaUnsmoothed = gradient(Mua);
D1Mua = sgolayfilt(D1MuaUnsmoothed, order, frameLen);
D2MuaUnsmoothed = gradient(D1Mua);
D2Mua = sgolayfilt(D2MuaUnsmoothed, order, frameLen);

WaterIndices = find(wavelengths >= 800 & wavelengths <= 850);
WaterCostFunctionD2 = @(c) norm(D2extinctions(WaterIndices,:) * c' - D2Mua(WaterIndices))^2;
[WaterConcD2, ~] = fminsearchbnd(WaterCostFunctionD2, ConcGuess, LowerBound, UpperBound);
WaterConcD2 = WaterConcD2(1);

HbIndices = find(wavelengths >= 753 & wavelengths <= 771);
HbCostFunctionD2 = @(c) norm(D2extinctions(HbIndices,:) * c' - D2Mua(HbIndices))^2;
[HbConcD2, ~] = fminsearchbnd(HbCostFunctionD2, ConcGuess, LowerBound, UpperBound);
HbD2Mua = D2extinctions * [HbConcD2(1), 0, 0]';
HbConcD2 = HbConcD2(2);



% DIRECT FITTING

    figure(4);
    h1 = plot(Mua, 'b');
    hold on;

    if TrueValue ~= 0
        h3 = plot(TrueValue, '--y');
        drawnow;
    end


%LowerBound(1:2) = [WaterConcD2(1), HbConcD2(2)];
%UpperBound(1:2) = [WaterConcD2(1), HbConcD2(2)];


lambda = 0;
costFunctionRegular = @(c) norm([extinctionsUnsmoothed * c' - AbsorbanceVector; lambda * (c' - ConcGuess')])^2;
costFunction = @(c) norm([extinctions * c' - Mua; lambda * (c' - ConcGuess')])^2;
D1costFunction = @(c) norm([D1extinctions * c' - D1Mua; lambda * (c' - ConcGuess')])^2;
D2costFunction = @(c) norm([D2extinctions * c' - D2Mua; lambda * (c' - ConcGuess')])^2;


Concentrations1 = fminsearchbnd(D2costFunction, ConcGuess, LowerBound, UpperBound);


ConcGuess2 = [WaterConcD2, ConcGuess(2), ConcGuess(3)];
LowerBound2 = [WaterConcD2, LowerBound(2), LowerBound(3)];
UpperBound2 = [WaterConcD2, UpperBound(2), UpperBound(3)];

    % Next, we solve for H2O. This is done with D1.

Concentrations2 = fminsearchbnd(D1costFunction, ConcGuess2, LowerBound2, UpperBound2);

ConcGuess3 = Concentrations2;

LowerBound3 = LowerBound;
UpperBound3 = UpperBound;

%LowerBound3 = [LowerBound2(1), Concentrations2(2:3)];
%UpperBound3 = [UpperBound2(1), Concentrations2(2:3)];

    % Finally, we solve for the rest.

Concentrations0 = fminsearchbnd(costFunctionRegular, ConcGuess3, LowerBound3, UpperBound3);

if DerivativeOrder == 0
    Concentrations = Concentrations0;
elseif DerivativeOrder == 1
    Concentrations = Concentrations2;
elseif DerivativeOrder == 2
    Concentrations = [Concentrations1(1), Concentrations2(2:3)];
elseif DerivativeOrder == -1
    Concentrations = [WaterConcD2, HbConcD2, Concentrations2(3)];
end

    h2 = plot(extinctions * Concentrations', 'r');
    hold off;
    drawnow;

else

AbsorbanceSmoothed = smoothdata(AbsorbanceVector, 1, "movmean", RegularFitting);


% DIRECT FITTING
    figure(4);
    h1 = plot(AbsorbanceSmoothed, 'b');
    hold on;

    if TrueValue ~= 0
        h3 = plot(TrueValue, '--y');
        drawnow;
    end

% Optimize the concentration estimates using a nonlinear solver

lambda = 0;
costFunction = @(c) norm([extinctionsUnsmoothed * c' - AbsorbanceVector; lambda * (c' - ConcGuess')])^2;

Concentrations = fminsearchbnd(costFunction, ConcGuess, LowerBound, UpperBound);

    h2 = plot(extinctionsUnsmoothed * Concentrations', 'r');
    hold off;
    drawnow;


end