function TPSF = TwoLayerTPSF(Mua1, Mua2, MuSprime1, MuSprime2, A, L, rho, times)

% This is equation 43 from Liemert and Kienle's 2010 paper for N-layered
% turbid media. It provides the Green's function, which in our case is the
% TPSF. z = 0 since we want reflectance.

MuSprime = MuSprime1;

t = times';

Ns = 5000;              % Number of points for the Integration. 
Nlambda = length(Mua1);
Ntime = length(times);

e1 = ones(Nlambda, 1);
e2 = ones(Ntime, 1);
e3 = ones(1, 1, Ns);

Reff=0.493;
D=1./(3 .* (MuSprime + (Mua1 + Mua2)/2));
z0=1./(Mua1+MuSprime);
zb = 2.*D.*((1 + Reff)/(1 - Reff));


c = 299.792/(1.4); % speed of light (mm/ns)

function IntegrandValue = IntegrandFFT(w)

    alpha1 = sqrt(Mua1 ./ D + 1i .* w ./ D ./ c); 
    alpha2 = sqrt(Mua2 ./ D + 1i .* w ./ D ./ c);
    Da = (alpha1 - alpha2) ./ (alpha1 + alpha2);

    IntegrandSection1 = cosh(alpha1.*(z0+2.*zb)) - cosh(alpha1.*z0); % Dimensions (Nlambda, 1)
    IntegrandSection2 = D .* alpha1 .* (Da + exp(2 .* alpha1 .* (L + zb))); % (Nlambda, 1)
    %ExpTerm = exp(1i .* (t .* w)); % (1, Ntime, Ns)

    IntegrandValue = (IntegrandSection1 .* Da) ./ (e1 .* IntegrandSection2); % 

    % ( (Nlambda, 1) x (1, Ntime) x ( (Nlambda, 1) x (1,Ntime) )
    % Total: (Nlambda, Ntime)
end


function IntegrandValue = Integrand(w)
    
    w = reshape(w, [1, 1, length(w)]);
    IntegrandSection3 = exp(1i .* (t' .* w)); % (1, Ntime, Ns)
    IntegrandValue = IntegrandFFT(w) .* IntegrandSection3;
end


%Nw = length(t)*2 - 2;
%wmax = 1/t(length(t));
%w = linspace(-wmax, wmax, Nw);
%dw = w(2) - w(1);

%Fw = IntegrandFFT(w);
%Fwshift = ifftshift(Fw);
%Ft = ifft(Fwshift, [], 2);

%dt = 2*pi / (Nw*dw);
%tvecfft = (-Nw/2:Nw/2-1) * dt;


%FtInterpolated = interp1(tvecfft(end/2+1:end), Ft(1,end/2+1:end), t, 'linear', 'extrap')';

%Integral = FtInterpolated;

wvec = linspace(-3e1, 3e1, Ns);

Integral = trapz(wvec, Integrand(wvec), 3); % Should be (Nlambda, Ntime, Ns) --> (Nlambda, Ntime)


    function Green = GreensFunction(t)

        InverseFourPiDct = (e1 * e2') ./ (4*pi*c*D*t'); % (Nlambda, 1) x (1, Ntime)
        GreenSection1 = c .* InverseFourPiDct.^(3/2) .* exp(-c.*Mua1*t') .* exp(-pi .* rho.^2 .* InverseFourPiDct); % (Nlambda, Ntime)
        GreenSection2 = exp(-((-z0).^2 * e2') .* pi.* InverseFourPiDct) - exp(-((z0+2.*zb).^2 * e2') .* pi .* InverseFourPiDct); % (Nlambda, 1) x (1, Ntime) x (Nlambda, Ntime)
        GreenSection3 = InverseFourPiDct .* exp( - rho.^2 .* pi .* InverseFourPiDct) .* (1 ./ (2.*pi));

        
        Green = (GreenSection1 .* GreenSection2) + (GreenSection3 .* Integral);

    end


    TPSF = A .* GreensFunction(t);
    TPSF = TPSF(:,1:end);

    %figure(27)
    %plot(times, TPSF)

end