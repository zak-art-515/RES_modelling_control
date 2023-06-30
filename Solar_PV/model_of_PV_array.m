%% Information from the KC200GT solar array datasheet
Iscn = 8.21; % Nominal short circuit voltage [A]
Vocn = 32.9; % Nominal array open circuit voltage [V]
Imp = 7.61; % Array current @ maximum power point [A]
Vmp = 26.3; % Array voltage @ maximum power point [V]
Pmax_e = Vmp*Imp; % Array maximum output peak power [W]
Kv = -0.123; % Voltage/temperature coefficient [V/K]
Ki = 3.18e-3; % Current/temperature coefficient [A/K]
Ns = 36; % Number of series cells

%% Constants
k = 1.3806503e-23; % Boltzmann [J/K]
q = 1.60217646e-19; % Electron Charge [Colomb]
a = 1.3; % Diode constant

%% Nominal values
Gn = 1000; % Nominal irradiance [W/m^2]
Tn = 25 + 273.15; % Nominal operating temperature [K]

%% Adjusting algorithm to nominal condition
G = 800;
T = 25 + 273.15;

Vtn = k * Tn / q; % Thermal junction voltage (nominal)
Vt = k * T / q; % Thermal junction voltage (current temperature)

Ion = Iscn/(exp(Vocn/a/Ns/Vtn)-1); % Nominal diode saturation current
Io = Ion;

% Reference values of Rs and Rp
Rs_max = (Vocn - Vmp)/Imp;
Rp_min = Vmp / (Iscn - Imp) - Rs_max;

%Initial values of Rp and Rs
Rp = Rp_min;
Rs = 0;

tol = 0.001; % Power mismatch Tolerance

P = (0);

error = Inf;

%% Iterative process for Rs and Rp until Pmax,model = Pmax,experimental

while error > tol
    % Temperature and irradiance effect on the current
    dT = T - Tn;
    Ipvn = (Rs + Rp)*Iscn/Rp; % nominal light-generated current
    Ipv = (Ipvn + Ki*dT)*G/Gn; % actual light-generated current
    Isc = (Iscn + Ki*dT)*G/Gn; % actual short-circuit current

    % Increments of Rs
    Rs = Rs + 0.01;

    % Parallel resistance
    Rp = Vmp*(Vmp + Imp*Rs)/(Vmp*Ipv - Vmp*Io*exp((Vmp + Imp*Rs)/Vt/Ns/a) + Vmp*Io - Pmax_e);

    % Solving the I-V equation for several (V,I) pairs
    clear V
    clear I

    V = 0:0.1:50; % Voltage vector
    I = zeros(1,size(V,2)); % current vector

    for j = 1:size(V,2) % compute for all voltage values
        % Solves g = I - f(I,V) = 0 with Newton Raphson method

        g(j) = Ipv - Io*(exp((V(j)+Rs*I(j))/Vt/Ns/a - 1)) -  (V(j) + Rs*I(j))/Rp - I(j);

        while abs(g(j)) > 0.001
            g(j) = Ipv - Io*(exp((V(j)+Rs*I(j))/Vt/Ns/a - 1)) -  (V(j) + Rs*I(j))/Rp - I(j);
            glin(j) = - Io*Rs/Vt/Ns/a*(exp((V(j)+Rs*I(j))/Vt/Ns/a)) -  Rs/Rp - 1;
            I_(j) = I(j) - g(j)/glin(j);
            I(j) = I_(j);
        end

    end

    % Calculates power using I-V equation
    P = (Ipv - Io*(exp((V + I.*Rs)/Vt/Ns/a)-1) - (V + I.*Rs)/Rp).*V;

    Pmax_m = max(P);

    error = (Pmax_m - Pmax_e);
end

%% Outpputs
fprintf('Model info:\n');
fprintf(' Rp_min = %f',Rp_min);
fprintf('\n Rp = %f',Rp);
fprintf('\n Rs_max = %f',Rs_max);
fprintf('\n Rs = %f',Rs);
fprintf('\n a = %f',a);
fprintf('\n T = %f',T - 273.15);
fprintf('\n G = %f',G);
fprintf('\n Pmax,m = %f (model)',Pmax_m);
fprintf('\n Pmax,e = %f (experimental)',Pmax_e);
fprintf('\n tol = %f',tol);
fprintf('\n error = %f',error);
fprintf('\n Ipv = %f',Ipv);
fprintf('\n Isc = %f',Isc);
fprintf('\n Ion = %f',Ion);

%% Plots of V-I curve and V-P curve
figure(1)
plot(V,I)
ylim([0 9])
figure(2)
plot(V,P)
ylim([0 200])
[Pmax,Indmax] = max(P)
V(Indmax)