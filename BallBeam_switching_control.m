% 시스템 파라미터
Kbb = 41.83;     % [cm/(rad*s^2)]
tau = 0.0248;    % [s]
Ki = 1.5286;     % [rad/V*s]
Lbeam = 42.55;   % [cm]
rarm = 2.54;     % [cm]
rb = 1.27;       % [cm]
mb = 0.064;      % [Kg]

% 자코비안 제어기 상태공간
Aj = [0 1 0 0; 0 0 Kbb 0; 0 0 0 1; 0 0 0 -1/tau];
Bj = [0; 0; 0; Ki / tau];
C = [1 0 0 0];

% 2단 선형화 상태공간
A2 = [0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 -1/tau];
B2 = [0; 0; 0; 1];

% 제어기/관측기 극점, 값은 같은데 이해하기 쉽게 구분
poles_K = [-1.7 -1.7 -1.7 -1.7];
poles_L = [-15 -25 -35 -45]';
Kj = acker(Aj, Bj, poles_K);
Lj = acker(Aj', C', poles_L)';
K2  = acker(A2,  B2,  poles_K);
L2  = acker(A2', C', poles_L)';

% 초기 상태
x0 = [-18; 0; -0.611; 0];
z1 = x0(1);
z2 = x0(2);
z3 = Kbb * sin(x0(3));
z4 = Kbb * cos(x0(3)) * x0(4);
zi = [z1; z2; z3; z4];
Z0 = [x0; zi];
tspan = [0 7];

% 시뮬레이션 실행
[T, Z] = ode45(@(t, Z) switching(t, Z, Kbb, tau, Ki, Kj, Lj, K2, L2), tspan, Z0);

% 결과 시각화
figure;
plot(T, Z(:, 1), 'k-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Ball Position [cm]');
title('Switching');
grid on;

% 실시간 스위칭 제어 함수
function dZ = switching(~, Z, Kbb, tau, Ki, Kj, Lj, K2, L2)
    x = Z(1:4);
    zhat = Z(5:8);
   
    y = x(1);

    if x(1) < -17.8
        % 자코비안 제어기
        xhat = zhat;
        yhat = zhat(1);
        u = -Kj * xhat;

        x3_hat = x(3);             
        x4_hat = x(4);
        A = [0 1 0 0; 0 0 Kbb 0; 0 0 0 1; 0 0 0 -1/tau];
        B = [0; 0; 0; Ki / tau];
        L = Lj;
        v = 0;
    else
        % 2단 선형화 제어기
        yhat = zhat(1);
        x3_hat = asin(zhat(3)/Kbb);
        x4_hat = zhat(4)/(Kbb * cos(x3_hat));
        xhat = [zhat(1); zhat(2); x3_hat; x4_hat];

        v = -K2 * zhat;
        u = (tau / (Ki * Kbb * cos(x3_hat))) * v;
        A = [0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 -1/tau];
        B = [0; 0; 0; 1];
        L = L2;
    end

    x1_dot = x(2);
    x2_dot = Kbb * sin(x3_hat);
    x3_dot = x4_hat;
    x4_dot = -x4_hat / tau + (Ki * u) / tau;
    dx = [x1_dot; x2_dot; x3_dot; x4_dot];
    dxhat = A * zhat + B * v + L * (y - yhat);

    dZ = [dx; dxhat];
end
