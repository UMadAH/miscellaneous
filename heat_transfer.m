% define parameters
Cl = 0.24;
h = 0.01;
k = 48;
rho = 7.86e+3;
Cv = 0.46e+3;
a = k / (rho * Cv);
lambda = Cl / a;
tau = lambda * h ^ 2;

T_inf = 298;
q_source = 50;
h_conv = 3;

% mesh
m = 100;
M = m + 1;
n = 200;
N = n + 1;

m0 = floor((M + 1)  / 2);
n0 = floor((N + 1) / 2);

T = ones(M, N) * T_inf;
T0 = repmat(T, 1);

% heat sources
i_type = 1;             % 0 as point
                             % 1 as linear
sigma = 0;


% iteration
max_time = 10000;
for i = 1:1:max_time
    for j = 2:1:(M - 1)
        for k = 2:1:(N - 1)
            delta_x = (T(j + 1, k) - 2 * T(j , k) + T(j - 1, k)) / (h * h);
            delta_y = (T(j, k + 1) - 2 * T(j , k) + T(j, k - 1)) / (h * h);
            if i_type == 0
                if (j == m0) && (k == n0)
                    T0(j, k) = T(j , k) + tau * a * (delta_x + delta_y) + tau * q_source / k;
                else
                    T0(j, k) = T(j , k) + tau * a * (delta_x + delta_y);
                end
            elseif i_type == 1
                T0(j, k) = T(j , k) + tau * a * (delta_x + delta_y) + tau * gauss(j, k, m0, n0, q_source) / k;
            end
        end
    end
    
    % L
    for k = 2:1:(N - 1)
        T0(1, k) = (h_conv * T_inf + (k / h) * T(2, k)) / (h_conv + (k / h));
    end
    
    % R
    for k = 2:1:(N - 1)
        T0(M, k) = (h_conv * T_inf - (k / h) * T(M - 1, k)) / (h_conv - (k / h));
    end
    
    % T
    for j = 1:1:(M)
        T0(j, 1) = (h_conv * T_inf + (k / h) * T(j, 2)) / (h_conv + (k / h));
    end
    
    % B
    for j = 1:1:(M)
        T0(j, N) = (h_conv * T_inf - (k / h) * T(j, N - 1)) / (h_conv - (k / h));
    end
    T = repmat(T0, 1);
    gca = pcolor(T0);
    set(gca, 'LineStyle','none');
    colorbar
    pause(0.01)
end


function y = gauss(x, y, x0, y0, q_source)
f_diff = (x - x0)^ 2 + (y - y0) ^ 2;
f_diff = sqrt(f_diff);
if f_diff <= 5
    y = (5 - f_diff) * q_source / 5;
else
    y = 0;
end
end