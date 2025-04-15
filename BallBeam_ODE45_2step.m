% 시스템 파라미터
Kbb = 41.83;     % [cm/(rad*s^2)]
tau = 0.0248;    % [s]
Ki = 1.5286;     % [rad/V*s]
Lbeam = 42.55;   % [cm]
rarm = 2.54;     % [cm]
rb = 1.27;       % [cm]
mb = 0.064;      % [Kg]

% 자코비안 선형화 모델
A0 = [0 1 0 0;
      0 0 1 0;
      0 0 0 1;
      0 0 0 -1/tau];

B0 = [0; 0; 0; 1];

C = [1 0 0 0];  % z1 = x1

% 제어기/관측기 극점
poles_K = [-1.7 -1.7 -1.7 -1.7]; % 제어기 -1.7 고정
poles_L = [-15 -25 -35 -45]'; % 관측기

% 이득 계산
K = acker(A0, B0, poles_K);
L = acker(A0', C', poles_L)';

% 초기 상태
x0 = [-18; 0; -0.611; 0];     % 실제 시스템 초기 상태

% 보조 상태 z 초기값 계산
z1 = x0(1);
z2 = x0(2);
z3 = Kbb * sin(x0(3));
z4 = Kbb * cos(x0(3)) * x0(4);
zi = [z1; z2; z3; z4];     % 관측기 초기값 (보조 상태 기준)

Z0 = [x0; zi];             % 8차 상태 초기화

% 시간축
tspan = [0 7];

% ODE45를 이용한 비선형 시스템 시뮬레이션
% 형식: [T, Z] = ode45(@(t, Z) 시스템함수(t, Z, 기타파라미터...), 시간범위, 초기값)
% ▷ T: 시간 벡터 (예: 0, 0.1, 0.2, ..., 7)
% ▷ Z: 각 시간에서의 상태값 (행 = 시간, 열 = 상태)
[T, Z] = ode45(@(t, Z) ballbeam2(t, Z, Kbb, tau, Ki, K, L), tspan, Z0);

% 결과 시각화
figure;
plot(T, Z(:, 1), 'r--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Ball position [cm]');
title('2단 선형화 기반 제어기 결과');
grid on;


function dZ = ballbeam2(~, Z, Kbb, tau, Ki, K, L)
    x = Z(1:4);      % 실제 상태
    zhat = Z(5:8);   % 추정 상태
    Z = [x;zhat];

    % 관측기 오차
    y = x(1);
    yhat = zhat(1);
    
    % z → x 역변환
    x3_hat = asin(zhat(3)/Kbb);
    x4_hat = zhat(4)/(Kbb * cos(x3_hat));
    xhat = [zhat(1); zhat(2); x3_hat; x4_hat];

    % 제어 입력 v (z 기반 제어기)
    v = -K * zhat;

    % 입력 변환 u (논문 식 27)
    u = (tau / (Ki * Kbb * cos(x3_hat))) * v;

    % 실제 시스템의 비선형 모델
    x1_dot = x(2);
    x2_dot = Kbb * sin(x3_hat);
    x3_dot = x4_hat;
    x4_dot = -x4_hat/tau + (Ki * u)/tau;
    dx = [x1_dot; x2_dot; x3_dot; x4_dot];

    % 자코비안 선형화 모델
    A0 = [0 1 0 0;
          0 0 1 0;
          0 0 0 1;
          0 0 0 -1/tau];
    B0 = [0; 0; 0; 1];
    C = [1 0 0 0];  

    % 관측기 (논문 식 23)
    dxhat = A0 * zhat + B0 * v + L * (y - yhat);
    dZ = [dx; dxhat];
end
