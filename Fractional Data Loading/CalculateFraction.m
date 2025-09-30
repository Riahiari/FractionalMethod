function f = CalculateFraction(AbsT, AbsRec, UseMean, PlotFraction)

%%%% f calculation %%%%%
if true
    f = zeros(size(AbsRec, 1), 1);
    for j=1:size(AbsRec, 1)
        MuaFit = polyfit(AbsT(j,:), AbsRec(j,:), 1);
        f(j) = abs(MuaFit(1));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%z


%%%%% PLOT FRACTION %%%%%
if PlotFraction
    figure(11)
    clf
    plot(f)
end
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% mean f %%%%%
if UseMean == 1

f = median(f);

elseif UseMean == 2

f = smoothdata(f, 1, "movmedian", 100);

elseif UseMean == 0

f = f;

end
%%%%%%%%%%%%%%%%%%%


end