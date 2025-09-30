function TPSF = OneLayerTPSF(mua, mus, A, rho, times)


% Calculate theoretical TPSF
    function Ref = ReflectanceOneLayer(t, mua, mus, A, rho)

        e1 = ones(length(mua), 1);
        e2 = ones(length(t), 1);

        v = 299.792; %speed of light in vaccum expressed in mm/ns
        n = 1.4; %refractive index of the medium
        c = v/n; %speed of light in the medium

        D = e1 ./ (3.*mus);%diffusion coefficient in mm
        zo = e1 ./ (mua+mus); %first point source term in mm
        zb = 2 * (1+0.493)/(1-0.493) .* D; %2nd point source term in mm
        

        alpha = -c .* (mua * t');
        beta = c * (4*pi*c .* D * t').^-1.5000;
        gamma = -(zo.^2 + e1 .* rho.^2) ./ (4*c .* D * t');
        delta = -((zo + 2*zb).^2 + e1 .* rho.^2) ./ (4*c .* D * t');
        epsilon = 1/2 * (4*pi*c .* D).^-1.5000 .* zo;
        phi = 1 + 2*zb ./ zo;
        
        Ref = A .* exp(alpha) .* (0.118 .* (beta .* (exp(gamma) - exp(delta)) ) + ...
              0.306 .* ( (epsilon * (t'.^-2.5000) ) .* (exp(gamma) + phi .* exp(delta))));

    end

    times = times';

    TPSF = ReflectanceOneLayer(times, mua, mus, A, rho); % Calculate TPSF using reflectance function

end