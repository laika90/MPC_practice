clear all;
clf;

z0 = 0;
z1 = 0;
x0  = [z0; z1];
x_size = length(x0);
x = x0;
ut = 0;

iter_num = 30;

result = zeros(x_size, iter_num);
time   = 1:iter_num;

for i = 1:iter_num
    result(:, i) = x;
    x = damped_oscilation(x, ut);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [1.750 -0.902;
     1     0;];
B = [1; 0];
C = [1  0];

% horizon
T = 10;

F = [];
G = [];
H = [];
current_Apower = eye(size(A));
zeromat_B = zeros(size(B));
zeromat_C = zeros(size(C));

for i=1:T
    % construct F
    current_Apower = current_Apower * A;
    F = [F; current_Apower];

    % construct G and H
    row_G = [];
    row_H = [];
    for j=1:T
        if i == j 
            row_G = [row_G, B];
            row_H = [row_H, C];
        elseif i > j
            row_G = [row_G, A^(i-j)*B];
            row_H = [row_H, zeromat_C];
        else 
            row_G = [row_G, zeromat_B];
            row_H = [row_H, zeromat_C];
        end
    end
    G = [G; row_G];
    H = [H; row_H];
end

% 多分同じ，修正予定
size_u = 1;
D = [];
eye_u = eye(size_u);
zeromat_u = zeros(size_u, size_u);

for i=1:T-1
    row_D = [];

    for j =1:T
        if i == j
            row_D = [row_D, -eye_u];
        elseif i+1 == j
            row_D = [row_D, eye_u];
        else
            row_D = [row_D, zeromat_u];
        end
    end

    D = [D; row_D];
end

% cost function coefficient
alpha_z = 0.1;
alpha_u = 0.5;
alpha_v = 1.0;

z_goal = 1;
y_ref = z_goal*ones(T, 1);

M = transpose(D) * D;
eye_M = eye(size(M));

Q = 2 * (transpose(G)*transpose(H)*(alpha_z*eye_M + alpha_v*M)*H*G + alpha_u*M);
% R  初期値xが反復ごとに変わるので，ループ内部で計算 


% constraints
Z = [];
for i=1:2*T
    row_Z = [];
    for j =1:T
        if mod(i, 2) == 1
            if fix((i+1)/2) == j
                row_Z = [row_Z, eye_u];
            else
                row_Z = [row_Z, zeromat_u];
            end
        else
            if fix(i/2) == j
                row_Z = [row_Z, -eye_u];
            else
                row_Z = [row_Z, zeromat_u];
            end
        end
    end
    Z = [Z; row_Z];
end


% ここらへんなんとかしたい
u_abs_max = 0.2;
w = u_abs_max*ones(2*T, 1);


result_mpc = zeros(x_size, iter_num);
result_U   = zeros(T, iter_num);
result_R = [];

x_mpc = x0;


for i = 1:iter_num
    result_mpc(:, i) = x_mpc;
    R = 2*transpose(G)*transpose(H)*((alpha_z*eye_M+alpha_v*transpose(M))*H*F*x_mpc - alpha_z*y_ref);
    [U,fval,exitflag,output,lambda] = quadprog(Q, R, Z, w);
    x_mpc = damped_oscilation(x_mpc, U(1));
    result_U(:, i) = U;
    result_R = [result_R, R];
end

hold on;

figure(1);
plot(time, result(1, :));
xlabel('time');
ylabel('x');
title('線形振動子');
grid on;

figure(2);
subplot(2,1,1);
plot(time, result_mpc(1, :));
xlabel('time');
ylabel('x');
title('振動子(MPC)');
grid on;

subplot(2,1,2);
plot(time, result_U(1, :));
xlabel('time');
ylabel('u');
title('制御量');
grid on;






