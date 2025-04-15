function ODE45

% 시스템 파라미터
Kbb = 41.83;     % [cm/(rad*s^2)]
tau = 0.0248;    % [s]
Ki = 1.5286;     % [rad/V*s]
Lbeam = 42.55;   % [cm]
rarm = 2.54;     % [cm]
rb = 1.27;       % [cm]
mb = 0.064;      % [Kg]

% 상태공간 자코비안 선형화
A0 = [0 1 0 0;
      0 0 Kbb 0;
      0 0 0 1;
      0 0 0 -1/tau];

B0 = [0; 0; 0; Ki / tau];

% 출력 행렬
C = [1 0 0 0];
   
% 제어기/관측기 극점
poles_K = [-1.7 -1.7 -1.7 -1.7]; % 제어기 -1.7 고정
poles_L = [-15 -25 -35 -45]'; % 관측기

% 이득 계산
K = acker(A0, B0, poles_K);
L = acker(A0', C', poles_L)';

% 초기 상태
x0 = [-18; 0; -0.611; 0];     % 실제 시스템 초기 상태
xhat0 = [0; 0; 0; 0];         % 추정 상태 초기값
Z0 = [x0; xhat0];             % 8차 상태 초기화

% 시간축
tspan = [0 7];

% ODE45를 이용한 비선형 시스템 시뮬레이션
% 형식: [T, Z] = ode45(@(t, Z) 시스템함수(t, Z, 기타파라미터...), 시간범위, 초기값)
% ▷ T: 시간 벡터 (예: 0, 0.1, 0.2, ..., 7)
% ▷ Z: 각 시간에서의 상태값 (행 = 시간, 열 = 상태)
[T, Z] = ode45(@(t, Z) ballbeam(t, Z, Kbb, tau, Ki, A0, B0, K, L), tspan, Z0);

% 결과 시각화
figure;
    plot(T, Z(:, 1), 'b-', 'LineWidth', 2);  % 볼 위치 x1
    xlabel('Time (s)');
    ylabel('Ball position [cm]');
    title('Ball and Beam System - ODE45');
    grid on;
end

function dZ = ballbeam(~, Z, Kbb, tau, Ki, A0, B0, K, L)
    x = Z(1:4);
    xhat = Z(5:8);

    u = -K * xhat;

    x1_dot = x(2);
    x2_dot = Kbb * sin(x(3));
    x3_dot = x(4);
    x4_dot = -x(4)/tau + (Ki * u)/tau;
    dx = [x1_dot; x2_dot; x3_dot; x4_dot];

    y = x(1);
    yhat = xhat(1);
    dxhat = A0 * xhat + B0 * u + L * (y - yhat);

    dZ = [dx; dxhat];
end