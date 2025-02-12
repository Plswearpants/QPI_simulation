function Isq = Isq_integral(s1, s2, b, n)
    % Integral over phi1 and phi2 with specified number of points n
    [phi1, phi2] = meshgrid(linspace(0, pi, n));
    integrand = cos(s1 * phi1) .* cos(s2 * phi2) ./ (b - cos(phi1) - cos(phi2));
    Isq = trapz(trapz(integrand)) * (pi / n)^2;
end
