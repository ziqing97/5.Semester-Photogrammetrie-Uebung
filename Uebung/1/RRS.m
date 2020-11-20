function[X0, Y0, Z0, omega, phi, kappa, Q] = RRS(c, X0, Y0, Z0, kappa, bildkoordinaten, koordinaten)
phi = 0;    % Horizontierung der Kamera durch stabilisierende Plattform 
omega = 0;
delta_x = zeros(6,1);
epsilon = 1;        % irgend eine grosse Zahl
x0 = 0;
y0 = 0;

i = 0;
while max(epsilon) >  1e-3
    X0 = X0 + delta_x(1);
    Y0 = Y0 + delta_x(2);
    Z0 = Z0 + delta_x(3);
    omega = omega + delta_x(4);
    phi = phi + delta_x(5);
    kappa = kappa + delta_x(6);

    r11 = cos(phi) * cos(kappa);
    r12 = -cos(phi) * sin(kappa);
    r13 = sin(phi);
    r21 = sin(omega) * sin(phi) * cos(kappa) + cos(omega) * sin(kappa);
    r22 = -sin(omega) * sin(phi) * sin(kappa) + cos(omega) * cos(kappa);
    r23 = -sin(omega) * cos(phi);
    r31 = -cos(omega) * sin(phi) * cos(kappa) + sin(omega) * sin(kappa);
    r32 = cos(omega) * sin(phi) * sin(kappa) + sin(omega) * cos(kappa);
    r33 = cos(omega) * cos(phi);

    N = r13 * (koordinaten(:,1) - X0) + r23 * (koordinaten(:,2) - Y0) + r33 * (koordinaten(:,3) - Z0);
    Zx = r11 * (koordinaten(:,1) - X0) + r21 * (koordinaten(:,2) - Y0) + r31 * (koordinaten(:,3) - Z0);
    Zy = r12 * (koordinaten(:,1) - X0) + r22 * (koordinaten(:,2) - Y0) + r32 * (koordinaten(:,3) - Z0);

    xn = x0 - c .* Zx ./ N;
    yn = y0 - c .* Zy ./ N;
    
    dxdX0 = -c ./ (N.^2) .* (r13 .* Zx - r11 .* N);
    dxdY0 = -c ./ (N.^2) .* (r23 .* Zx - r21 .* N);
    dxdZ0 = -c ./ (N.^2) .* (r33 .* Zx - r31 .* N);
    dxdomega = -c ./ N .* (((koordinaten(:,2) - Y0) .* r33 -    ... 
              (koordinaten(:,3) - Z0) .* r23) .* Zx ./ N - (koordinaten(:,2) - Y0) .* r31 + (koordinaten(:,3) - Z0) .* r21);
    dxdphi = c ./ N .* ((Zx .* cos(kappa) - Zy .* sin(kappa)) .* Zx ./ N + N .* cos(kappa));
    dxdkappa = -c ./ N .* Zy;

    dydX0 = -c ./ (N.^2) .* (r13 .* Zy - r12 .* N);
    dydY0 = -c ./ (N.^2) .* (r23 .* Zy - r22 .* N);
    dydZ0 = -c ./ (N.^2) .* (r33 .* Zy - r32 .* N);
    dydomega = -c ./ N .* (((koordinaten(:,2) - Y0) .* r33 -    ... 
              (koordinaten(:,3) - Z0) .* r23) .* Zy ./ N - (koordinaten(:,2) - Y0) .* r32 + (koordinaten(:,3) - Z0) .* r22);
    dydphi = c ./ N .* ((Zx .* cos(kappa) - Zy .* sin(kappa)) .* Zy ./ N - N .* sin(kappa));
    dydkappa = c ./ N .* Zx;

    A = [dxdX0, dxdY0, dxdZ0, dxdomega, dxdphi, dxdkappa;
         dydX0, dydY0, dydZ0, dydomega, dydphi, dydkappa];
    delta_x = (A' * A) \ A' * [bildkoordinaten(:,1) - xn; bildkoordinaten(:,2) - yn];

    epsilon = [abs(delta_x(1)), abs(delta_x(2)), abs(delta_x(3))];
    i = i + 1;
end
    y_dach = A * delta_x;
    e_dach = [bildkoordinaten(:,1) - xn; bildkoordinaten(:,2) - yn] - y_dach;
    sigma = (e_dach' * e_dach / (length(y_dach) - 6));
    Q = inv(A' * A);   % Inversion der Normalgleichungsmatrix, Gewicht;
    Q = sigma * Q;
end