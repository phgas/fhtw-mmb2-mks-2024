function simulate_system()
    G = 9.81;                            % Gravity acceleration [m/s^2]
    PHI_0_DEG = 30;                      % Angle of tilted carousel [°]
    PHI_0_RAD = deg2rad(PHI_0_DEG);      % Angle of tilted carousel [radians]
    R_KARUSELL = 6;                      % Radius of carousel [m]
    C = 20000;                           % Spring stiffness [N/m]
    M_GONDEL = 300;                      % Mass of gondola [kg]
    N_GONDEL = 0.2333;                   % Rotational speed of gondola [radians/second]
    B_GONDEL = 1.5;                      % Width of gondola [m]
    H_GONDEL = 1.5;                      % Height of gondola [m]
    DELTA_L = 0.005;                     % Max length of locking mechanism [m]
    
    % Initial conditions [x, x_dot, alpha, alpha_dot]
    x0 = R_KARUSELL;
    x_dot0 = 0;
    alpha0 = 0;
    alpha_dot0 = 2 * pi * N_GONDEL;
    
    % Time span
    tspan = [0, 30];
    
    % ODE solver
    [t, y] = ode45(@(t, y) odefun(t, y, M_GONDEL, C, G, PHI_0_RAD, DELTA_L, B_GONDEL, H_GONDEL, R_KARUSELL), tspan, [x0, x_dot0, alpha0, alpha_dot0]);
    
    figure;
    subplot(4, 1, 1);
    plot(t, y(:, 1), 'DisplayName', '$x(t)$', 'LineWidth', 1.5);
    ylabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 14);
    legend('Interpreter', 'latex', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 14);
    
    subplot(4, 1, 2);
    plot(t, y(:, 2), 'DisplayName', '$\dot{x}(t)$', 'LineWidth', 1.5);
    ylabel('$\dot{x}$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14);
    legend('Interpreter', 'latex', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 14);
    
    subplot(4, 1, 3);
    plot(t, y(:, 3), 'DisplayName', '$\alpha(t)$', 'LineWidth', 1.5);
    ylabel('$\alpha$ [rad]', 'Interpreter', 'latex', 'FontSize', 14);
    legend('Interpreter', 'latex', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 14);
    
    subplot(4, 1, 4);
    plot(t, y(:, 4), 'DisplayName', '$\dot{\alpha}(t)$', 'LineWidth', 1.5);
    ylabel('$\dot{\alpha}$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 14);
    legend('Interpreter', 'latex', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 14);

    han = axes(gcf, 'visible', 'off');
    han.XLabel.Visible = 'on';
    xlabel(han, 'Time (t)', 'FontSize', 14);

    % Get local extremes
    get_local_extremes(y(:, 1), 'x');
    get_local_extremes(y(:, 2), 'x\_dot');
    get_local_extremes(y(:, 3), 'alpha');
    get_local_extremes(y(:, 4), 'alpha\_dot');
end

function dydt = odefun(~, y, M_GONDEL, C, G, PHI_0_RAD, DELTA_L, B_GONDEL, H_GONDEL, R_KARUSELL)
    x = y(1);
    x_dot = y(2);
    alpha = y(3);
    alpha_dot = y(4);
    
    x_ddot = x * alpha_dot^2 - (C / M_GONDEL) * (x + DELTA_L) + G * sin(PHI_0_RAD) * cos(alpha);
    alpha_ddot = -(2 * x_dot * alpha_dot + G * sin(PHI_0_RAD) * sin(alpha) * x) / (x^2 + (5 / 3) * (B_GONDEL^2 + H_GONDEL^2) + 20 * R_KARUSELL^2);
    
    dydt = [x_dot; x_ddot; alpha_dot; alpha_ddot];
end

function get_local_extremes(data_series, name)
    local_maxima = [];
    local_minima = [];

    for i = 2:length(data_series)-1
        if data_series(i) > data_series(i-1) && data_series(i) > data_series(i+1)
            local_maxima = [local_maxima, data_series(i)];
        end
        if data_series(i) < data_series(i-1) && data_series(i) < data_series(i+1)
            local_minima = [local_minima, data_series(i)];
        end
    end

    fprintf('Local maxima of %s: %s\n', name, mat2str(local_maxima));
    fprintf('Local minima of %s: %s\n', name, mat2str(local_minima));
end
