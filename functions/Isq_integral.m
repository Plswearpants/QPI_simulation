function Isq = Isq_integral(s1, s2, b, n)
    % Computes the double integral for the Green's function calculation
    %
    % Inputs:
    %   s1 - x-component of the normalized separation vector (s_x/a)
    %        Can be a scalar or array (e.g., size [M,N] for a grid of points)
    %   s2 - y-component of the normalized separation vector (s_y/a)
    %        Must have the same dimensions as s1 (e.g., size [M,N])
    %   b  - normalized energy parameter: (omega + i*epsilon - E0)/(2*t), scalar
    %   n  - number of discretization points for numerical integration, scalar integer
    %
    % Output:
    %   Isq - result of the double integral
    %         Has the same dimensions as s1 and s2 (e.g., size [M,N])
    %         For example, if s1 and s2 are [100,100] arrays, Isq will be [100,100]
    %         ∫∫ cos(s1*phi1)*cos(s2*phi2)/(b + cos(phi1) + cos(phi2)) dphi1 dphi2
    %         integrated over [0,π]×[0,π]
    
    [phi1, phi2] = meshgrid(linspace(0, pi, n));
    integrand = cos(s1 * phi1) .* cos(s2 * phi2) ./ (b + cos(phi1) + cos(phi2));
    Isq = trapz(trapz(integrand)) * (pi / n)^2;
end
