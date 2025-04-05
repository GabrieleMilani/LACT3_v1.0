function [polyL, Int_i] = Material(polyL, Int_i, c, phi, s, rho)
    % rho=rho*100*10^-9; % kg/mm3
    rho=rho*10^-9; % kg/mm3
    tanphi=tan(phi*pi/180);

    polyL(:, size(polyL,2)-1) = s;
    polyL(:, size(polyL,2)) = rho;

    Int_i(:,12) = c;
    Int_i(:,13) = tanphi;
end