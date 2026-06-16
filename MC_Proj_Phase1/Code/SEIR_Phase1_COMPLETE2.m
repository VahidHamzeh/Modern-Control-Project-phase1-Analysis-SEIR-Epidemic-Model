
clc; clear; close all;

%% =========================================================
%  SEC 0: USER SETTINGS  (edit these two lines only)
% =========================================================
SAVE_FIGS  = false;                       % true = save PNGs to disk
OUTPUT_DIR = 'D:\phase 1 modern\Figures'; % folder must exist (created auto)

if SAVE_FIGS
    if ~exist(OUTPUT_DIR, 'dir')
        mkdir(OUTPUT_DIR);
        fprintf('Created output folder: %s\n', OUTPUT_DIR);
    end
end

fprintf('=================================================================\n');
fprintf('  SEIR EPIDEMIC MODEL — PHASE 1 COMPLETE ANALYSIS\n');
fprintf('  Student: Amirreza Farzaneh  |  ID: 401102231\n');
fprintf('  Sharif University of Technology — Dr. Pourshamsi\n');
fprintf('  Reference: Soulaimani & Kaddar, IEEE Access 2023\n');
fprintf('=================================================================\n\n');

%% =========================================================
%  SEC 1: TWO INCIDENCE FUNCTIONS  (g1 and g2)
% =========================================================
% g1: Saturated incidence (Article reference) -- captures crowding effect
% g2: Bilinear incidence  (classical comparison baseline)
%
% Both give identical linearized A matrix at DFE (because I=0 removes difference)

g1 = @(S, I, beta, c) beta .* S .* I ./ (1 + c .* I.^2);   % Eq. 6
g2 = @(S, I, beta, ~) beta .* S .* I;                        % Eq. 7

fprintf('[SEC 1] Incidence functions defined:\n');
fprintf('  g1(S,I) = beta*S*I / (1 + c*I^2)   [Saturated -- Article]\n');
fprintf('  g2(S,I) = beta*S*I                  [Bilinear  -- Comparison]\n\n');

%% =========================================================
%  SEC 2: TWO PARAMETER SETS  (Tables 1 & 2 from article)
% =========================================================
% Set 1 (v=0.1): WITH vaccination -> intended R0 <= 1 scenario
% Set 2 (v=0)  : WITHOUT vaccination -> R0 > 1 endemic scenario
% Control case (v=4.0): our stable operating point for RMSE analysis

% --- Common parameters ---
mu      = 0.003;   % Natural death rate (1/day)
c_par   = 4;       % Incidence saturation coefficient
sigma   = 0.004;   % Exposed -> Infected rate (1/day)
delta   = 0.001;   % Disease-induced death rate, Exposed (1/day)
epsilon = 0.03;    % Immunity rate, Exposed (1/day)
gamma   = 0.07;    % Treatment/recovery rate (1/day)
d_inf   = 0;       % Disease-induced death rate, Infected (1/day)
alpha   = 0.94;    % Fractional order (reference only -- we use alpha=1)

% --- Set 1: v=0.1, beta=0.01, b=40 ---
b1   = 40;
beta1 = 0.01;
v1   = 0.1;

% --- Set 2: v=0, beta=0.9, b=30 ---
b2    = 30;
beta2 = 0.9;
v2    = 0;

% --- Control case: v=4.0, beta=0.01, b=40 (gives R0<1 for valid RMSE) ---
b_c   = 40;
beta_c = 0.01;
v_c   = 4.0;

fprintf('[SEC 2] Parameter sets loaded:\n');
fprintf('  Set 1 : b=%d, beta=%.3f, v=%.1f\n', b1, beta1, v1);
fprintf('  Set 2 : b=%d, beta=%.1f, v=%.1f\n', b2, beta2, v2);
fprintf('  Control: b=%d, beta=%.3f, v=%.1f\n\n', b_c, beta_c, v_c);

%% =========================================================
%  SEC 3: EQUILIBRIA AND BASIC REPRODUCTION NUMBER R0
% =========================================================
% Disease-Free Equilibrium: E0 = (S0, 0, 0, 0)
% S0 = b / (mu + v)

S0_1 = b1 / (mu + v1);     % Set 1
S0_2 = b2 / (mu + v2);     % Set 2
S0_c = b_c / (mu + v_c);   % Control

% R0 via next-generation matrix (van den Driessche & Watmough 2002)
% R0 = sigma * (dg/dI)|_E0 / [(sigma+mu+delta+epsilon)*(mu+gamma+d)]
%    = sigma * beta * S0 / [...]

denom_R0 = (sigma + mu + delta + epsilon) * (mu + gamma + d_inf);

R0_1 = (sigma * beta1 * S0_1) / denom_R0;
R0_2 = (sigma * beta2 * S0_2) / denom_R0;
R0_c = (sigma * beta_c * S0_c) / denom_R0;

fprintf('[SEC 3] Disease-Free Equilibria and R0:\n');
fprintf('  Set 1 : S0 = %8.4f  |  R0 = %.4f', S0_1, R0_1);
if R0_1 <= 1, fprintf('  --> Stable DFE\n'); else, fprintf('  --> Unstable DFE\n'); end
fprintf('  Set 2 : S0 = %8.2f  |  R0 = %.2f', S0_2, R0_2);
if R0_2 <= 1, fprintf('  --> Stable DFE\n'); else, fprintf('  --> Unstable DFE\n'); end
fprintf('  Control: S0 = %8.4f  |  R0 = %.4f', S0_c, R0_c);
if R0_c <= 1, fprintf('  --> Stable DFE\n\n'); else, fprintf('  --> Unstable DFE\n\n'); end

% --- Endemic equilibrium E** (exists only when R0 > 1) ---
% I** found by solving K(I) = 0 numerically
% K(I) = g(S**(I), I) / I - (sigma+mu+delta+epsilon)*(mu+gamma+d) / sigma

fprintf('[SEC 3b] Endemic Equilibrium E** for Set 2 (R0 >> 1):\n');
K_fn = @(I) compute_K(I, b2, mu, v2, beta2, c_par, sigma, delta, epsilon, gamma, d_inf);

% Bracket for root search
I_upper = b2 * sigma / ((sigma+mu+delta+epsilon)*(mu+gamma+d_inf)) - 1e-6;
try
    I_star = fzero(K_fn, [1e-6, I_upper]);
    E_star = (mu + gamma + d_inf) / sigma * I_star;
    S_star = b2/(mu+v2) - (sigma+mu+delta+epsilon)*(mu+gamma+d_inf)/(sigma*(mu+v2)) * I_star;
    R_star = (v2*S_star + gamma*I_star + epsilon*E_star) / mu;
    fprintf('  S** = %8.4f\n', S_star);
    fprintf('  E** = %8.4f\n', E_star);
    fprintf('  I** = %8.4f\n', I_star);
    fprintf('  R** = %8.4f\n\n', R_star);
catch
    fprintf('  (fzero failed -- Set 2 dynamics diverge rapidly)\n\n');
    S_star = NaN; E_star = NaN; I_star = NaN; R_star = NaN;
end

%% =========================================================
%  SEC 4: LINEARIZATION  ->  A, B, C, D
% =========================================================
% Symbolic form at DFE (Eq. 24):
%   A = [-(mu+v),  0,          -beta*S0,         0   ]
%       [0,        -(eps+sig+mu+del), +beta*S0,   0   ]
%       [0,         sigma,      -(mu+gamma+d),    0   ]
%       [v,         epsilon,     gamma,           -mu  ]
%
% Note: Block lower-triangular structure -> one eigenvalue ALWAYS = -mu
%
% Numerical Jacobian (central difference) -- works for BOTH g1 and g2
% dg/dS|_E0 = 0  (since I=0)
% dg/dI|_E0 = beta * S0  (same for g1 and g2 at I=0)

fprintf('[SEC 4] Linearization at DFE:\n');

% Build A,B,C,D for each scenario
[A1, B1, C1, D1] = linearize_SEIR([S0_1;0;0;0], beta1, v1, mu, epsilon, sigma, delta, gamma, d_inf);
[A2, B2, C2, D2] = linearize_SEIR([S0_2;0;0;0], beta2, v2, mu, epsilon, sigma, delta, gamma, d_inf);
[Ac, Bc, Cc, Dc] = linearize_SEIR([S0_c;0;0;0], beta_c, v_c, mu, epsilon, sigma, delta, gamma, d_inf);

n = size(Ac, 1);   % system order = 4

fprintf('\n  A (Control, v=4.0):\n');
disp(Ac);
fprintf('  B (Control): [%.4f; 0; 0; %.4f]\n', Bc(1), Bc(4));
fprintf('  C = [0 0 1 0]   D = 0\n\n');

% Build state-space objects
sys1 = ss(A1, B1, C1, D1);
sys2 = ss(A2, B2, C2, D2);
sys_c = ss(Ac, Bc, Cc, Dc);

%% =========================================================
%  SEC 5: EIGENVALUE ANALYSIS AND STABILITY
% =========================================================
eigs1 = sort(real(eig(A1)));
eigs2 = sort(real(eig(A2)));
eigs_c = sort(real(eig(Ac)));

fprintf('[SEC 5] Eigenvalues:\n');
fprintf('  Set 1  (v=0.1): [%.5f, %.5f, %.5f, %.5f]', eigs1(1),eigs1(2),eigs1(3),eigs1(4));
if all(eigs1 < 0), fprintf('  STABLE\n'); else, fprintf('  UNSTABLE\n'); end

fprintf('  Set 2  (v=0  ): [%.5f, %.5f, %.5f, %.5f]', eigs2(1),eigs2(2),eigs2(3),eigs2(4));
if all(eigs2 < 0), fprintf('  STABLE\n'); else, fprintf('  UNSTABLE\n'); end

fprintf('  Control (v=4.0): [%.5f, %.5f, %.5f, %.5f]', eigs_c(1),eigs_c(2),eigs_c(3),eigs_c(4));
if all(eigs_c < 0), fprintf('  STABLE\n\n'); else, fprintf('  UNSTABLE\n\n'); end

%% =========================================================
%  SEC 6: TRANSFER FUNCTION, POLES, ZEROS
% =========================================================
[num_tf, den_tf] = ss2tf(Ac, Bc, Cc, Dc);
tf_obj = tf(num_tf, den_tf);
sys_poles_c = pole(sys_c);
sys_zeros_c = tzero(sys_c);

fprintf('[SEC 6] Transfer Function G(s) = C(sI-A)^{-1}B  [Control case]:\n');
disp(tf_obj);
fprintf('  Poles: '); disp(sys_poles_c');
fprintf('  Zeros: '); disp(sys_zeros_c');

%% =========================================================
%  SEC 7: PARAMETRIC SENSITIVITY  (beta +/-5%, +/-10%)
% =========================================================
fprintf('[SEC 7] Sensitivity Analysis -- beta +/-5%%, +/-10%%:\n');
fprintf('  %-6s  %-8s  %-12s  %-12s  %-12s  %-12s  %-8s  %-10s\n',...
        'Change','beta','lambda1','lambda2','lambda3','lambda4','R0','Stable?');

beta_mults = [0.90, 0.95, 1.00, 1.05, 1.10];
pct_labels = {'-10%','-5%','0%','+5%','+10%'};
all_eigs_sens = zeros(4, length(beta_mults));
R0_sens       = zeros(1, length(beta_mults));

for k = 1:length(beta_mults)
    bv   = beta_c * beta_mults(k);
    dg_v = bv * S0_c;
    Av   = [-(mu+v_c),  0,                          -dg_v,           0;
             0,         -(epsilon+sigma+mu+delta),   dg_v,            0;
             0,          sigma,                     -(mu+gamma+d_inf), 0;
             v_c,        epsilon,                    gamma,           -mu];
    ev      = sort(real(eig(Av)));
    all_eigs_sens(:,k) = ev;
    R0_sens(k) = (sigma * dg_v) / denom_R0;
    stab_str = 'YES';
    if any(ev >= 0), stab_str = '*** NO ***'; end
    fprintf('  %-6s  %-8.5f  %-12.6f  %-12.6f  %-12.6f  %-12.6f  %-8.4f  %-10s\n',...
            pct_labels{k}, bv, ev(1), ev(2), ev(3), ev(4), R0_sens(k), stab_str);
end
fprintf('\n');

%% =========================================================
%  SEC 8: CANONICAL FORMS
% =========================================================
% 8a. Characteristic polynomial
char_poly = poly(Ac);
a3 = char_poly(2);  % s^3 coefficient
a2 = char_poly(3);  % s^2 coefficient
a1_cp = char_poly(4);  % s^1 coefficient
a0 = char_poly(5);  % s^0 (constant)

fprintf('[SEC 8] Characteristic Polynomial:\n');
fprintf('  p(s) = s^4 + (%.6f)s^3 + (%.6f)s^2 + (%.6f)s + (%.8f)\n',...
        a3, a2, a1_cp, a0);
fprintf('  a3=%.6f, a2=%.6f, a1=%.6f, a0=%.8f\n\n', a3, a2, a1_cp, a0);

% 8b. Jordan Canonical Form
[V_jord, J_jord] = eig(Ac);       % eig gives diagonal for distinct eigenvalues
B_jordan = V_jord \ Bc;            % Bj = V^{-1} * B
C_jordan = Cc * V_jord;            % Cj = C * V
fprintf('[SEC 8b] Jordan Form (diagonal -- all eigenvalues distinct):\n');
fprintf('  J = diag(%.6f, %.6f, %.6f, %.6f)\n',...
        J_jord(1,1), J_jord(2,2), J_jord(3,3), J_jord(4,4));
fprintf('  Bj = V^{-1}*B:  [%.4f; %.4f; %.4f; %.4f]\n',...
        B_jordan(1), B_jordan(2), B_jordan(3), B_jordan(4));
fprintf('  Cj = C*V:       [%.4f  %.4f  %.4f  %.4f]\n\n',...
        C_jordan(1), C_jordan(2), C_jordan(3), C_jordan(4));

% 8c. Controller Canonical Form (CCF)
% Standard form (Ogata/Chen):
%   Ac_ccf last row = [-a0, -a1, -a2, -a3]
%   Bc_ccf = [0; 0; 0; 1]   <-- STANDARD CCF input vector
%   Cc_ccf = [b0, b1, b2, b3]  (TF numerator coefficients)
Ac_ccf = [0,   1,    0,     0;
           0,   0,    1,     0;
           0,   0,    0,     1;
          -a0, -a1_cp, -a2, -a3];

Bc_ccf = [0; 0; 0; 1];              % STANDARD CCF Bc vector
Cc_ccf = [1, 0, 0, 0];             % CCF output (reading first state)

fprintf('[SEC 8c] Controller Canonical Form (CCF):\n');
fprintf('  Ac_ccf last row: [%.8f, %.6f, %.6f, %.6f]\n', -a0, -a1_cp, -a2, -a3);
fprintf('  Bc_ccf = [0; 0; 0; 1]  (standard CCF input vector)\n');
fprintf('  Ac_ccf =\n');
disp(Ac_ccf);

% 8d. Observer Canonical Form (OCF)
% Standard form (Ogata/Chen):
%   Ao_ocf last COLUMN = [-a0; -a1; -a2; -a3]
%   Co_ocf = [0, 0, 0, 1]   <-- STANDARD OCF output vector
%   Bo_ocf = [1; 0; 0; 0]
Ao_ocf = [0,  0,  0,  -a0;
           1,  0,  0,  -a1_cp;
           0,  1,  0,  -a2;
           0,  0,  1,  -a3];

Bo_ocf = [1; 0; 0; 0];
Co_ocf = [0, 0, 0, 1];              % STANDARD OCF output vector

fprintf('[SEC 8d] Observer Canonical Form (OCF):\n');
fprintf('  Ao_ocf last column: [-%.8f; -%.6f; -%.6f; -%.6f]\n', a0, a1_cp, a2, a3);
fprintf('  Co_ocf = [0, 0, 0, 1]  (standard OCF output vector)\n');
fprintf('  Ao_ocf =\n');
disp(Ao_ocf);

% Verify CCF and OCF have correct eigenvalues
eigs_ccf = sort(real(eig(Ac_ccf)));
eigs_ocf = sort(real(eig(Ao_ocf)));
fprintf('  CCF eigenvalues match original? %s\n', ...
        mat2str(all(abs(eigs_ccf - eigs_c) < 1e-8)));
fprintf('  OCF eigenvalues match original? %s\n\n', ...
        mat2str(all(abs(eigs_ocf - eigs_c) < 1e-8)));

%% =========================================================
%  SEC 9: STRUCTURAL ANALYSIS -- 5 METHODS
% =========================================================
fprintf('[SEC 9] Structural Analysis (5 methods, Control case v=4.0):\n\n');

% ---------- 9.1 METHOD 1: Controllability / Observability Matrices ----------
Co_mat = ctrb(Ac, Bc);
Ob_mat = obsv(Ac, Cc);
rank_Co = rank(Co_mat);
rank_Ob = rank(Ob_mat);

fprintf('  Method 1 -- Matrix rank:\n');
fprintf('    rank(Co) = %d / %d  -->  ', rank_Co, n);
if rank_Co == n, fprintf('FULLY CONTROLLABLE\n');
else, fprintf('PARTIALLY CONTROLLABLE\n'); end
fprintf('    rank(Ob) = %d / %d  -->  ', rank_Ob, n);
if rank_Ob == n, fprintf('FULLY OBSERVABLE\n\n');
else, fprintf('PARTIALLY OBSERVABLE  (y=I only)\n\n'); end

% With all 4 outputs
Ob_full = obsv(Ac, eye(4));
fprintf('    rank(Ob) with y=[S,E,I,R] = %d / %d  -->  ', rank(Ob_full), n);
if rank(Ob_full)==n, fprintf('FULLY OBSERVABLE\n\n'); end

% ---------- 9.2 METHOD 2: PBH Test ----------
fprintf('  Method 2 -- PBH Test (Popov-Belevitch-Hautus):\n');
fprintf('    Controllable <=> rank[A - lambda*I | B] = n  for all eigenvalues\n');
eigs_all = eig(Ac);
for k = 1:n
    lam = eigs_all(k);
    M_ctrl = [Ac - lam*eye(n),  Bc];
    M_obs  = [Ac - lam*eye(n); Cc];
    r_ctrl = rank(M_ctrl);
    r_obs  = rank(M_obs);
    ctrl_ok = ''; if r_ctrl < n, ctrl_ok = ' <-- FAILS'; end
    obs_ok  = ''; if r_obs  < n, obs_ok  = ' <-- FAILS'; end
    fprintf('    lambda=%.5f:  rank[A-lI|B]=%d%s  |  rank[A-lI;C]=%d%s\n',...
            real(lam), r_ctrl, ctrl_ok, r_obs, obs_ok);
end
fprintf('\n');

% ---------- 9.3 METHOD 3: Jordan Form ----------
fprintf('  Method 3 -- Jordan Form:\n');
fprintf('    Controllable <=> no row of Bj=V^{-1}B (for Jordan 1st rows) is zero\n');
B_j_check = abs(B_jordan);
for k = 1:n
    if B_j_check(k) < 1e-10
        fprintf('    Bj(%d) = %.2e  --> Mode %d NOT controllable\n', k, B_j_check(k), k);
    else
        fprintf('    Bj(%d) = %.4f  --> Mode %d controllable\n', k, B_j_check(k), k);
    end
end
fprintf('\n');

% ---------- 9.4 METHOD 4: Gramian ----------
fprintf('  Method 4 -- Controllability & Observability Gramians:\n');
fprintf('    A*Wc_gram + Wc_gram*A'' + B*B'' = 0  (Lyapunov equation)\n');
try
    Wc_gram = lyap(Ac, Bc*Bc');
    Wo_gram = lyap(Ac', Cc'*Cc);
    rank_Wc_gram = rank(Wc_gram);
    rank_Wo_gram = rank(Wo_gram);
    fprintf('    rank(Wc_gram) = %d / %d  -->  ', rank_Wc_gram, n);
    if rank_Wc_gram==n, fprintf('Fully Controllable\n');
    else, fprintf('Partially Controllable\n'); end
    fprintf('    rank(Wo_gram) = %d / %d  -->  ', rank_Wo_gram, n);
    if rank_Wo_gram==n, fprintf('Fully Observable\n\n');
    else, fprintf('Partially Observable\n\n'); end
catch ME
    fprintf('    lyap() requires stable system. Error: %s\n\n', ME.message);
    Wc_gram = NaN; Wo_gram = NaN;
end

% ---------- 9.5 METHOD 5: SVD ----------
fprintf('  Method 5 -- Singular Value Decomposition (SVD):\n');
sv_c = svd(Co_mat);
sv_o = svd(Ob_mat);
fprintf('    Singular values of Co: ');
fprintf('%.4f  ', sv_c'); fprintf('\n');
fprintf('    Singular values of Ob: ');
fprintf('%.4f  ', sv_o'); fprintf('\n');
fprintf('    Controllable <=> sigma_min(Co) > 0 : %s\n', mat2str(min(sv_c) > 1e-10));
fprintf('    Observable   <=> sigma_min(Ob) > 0 : %s\n\n', mat2str(min(sv_o) > 1e-10));

%% =========================================================
%  SEC 10: MINIMAL REALIZATION
% =========================================================
sys_min = minreal(sys_c);
n_min   = size(sys_min.A, 1);
fprintf('[SEC 10] Minimal Realization:\n');
fprintf('  Original order  : %d\n', n);
fprintf('  Minimal order   : %d\n', n_min);
if n_min == n
    fprintf('  --> System is already minimal\n\n');
else
    fprintf('  --> Reduced from order %d to %d\n\n', n, n_min);
end

%% =========================================================
%  SEC 11: SIMULATION  (Nonlinear vs Linearized)
%
%  TWO CONSISTENT SCENARIOS -- each uses the SAME operating point
%  for both nonlinear (ode45) and linear (lsim):
%
%  SCENARIO A (v=4.0, Control, RMSE validation):
%    - Linear system sys_c built at v0=4.0, S0=9.99
%    - x0 is small perturbation around THAT equilibrium
%    - Step deviation: Dv=+0.2 around v=4.0
%    - Sinusoidal deviation: 0.15*sin(0.1*t) around v=4.0
%
%  SCENARIO B (v=0.1, Groupmate's exact parameters):
%    - Linear system sys1 built at v0=0.1, S0=388.35
%    - x0 = [0.95*S0_1, 2, 5, 0] (groupmate's exact values)
%    - Step: v 0.1 -> 0.2  (deviation +0.1 around v=0.1)
%    - Sinusoidal: v(t)=0.1+0.05*sin(0.02t) (deviation around v=0.1)
%    - Note: DFE is unstable at v=0.1 (R0=5.6>1), so RMSE is larger
% =========================================================
fprintf('[SEC 11] Running simulations...\n\n');

t_end  = 200;
t_eval = linspace(0, t_end, 5000)';
mask   = t_eval > 30;

ode_opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

% ============================================================
%  SCENARIO A: Control case (v=4.0) -- for RMSE validation
% ============================================================
fprintf('  --- Scenario A: Control case (v=4.0, R0=0.1441 < 1) ---\n');

% Correct initial conditions: small perturbation around CONTROL equilibrium
R0_eq_c  = v_c * S0_c / mu;               % R equilibrium at v=4.0
x0_A     = [S0_c*0.995; 2; 0.5; R0_eq_c*0.995];
dx0_A    = x0_A - [S0_c; 0; 0; R0_eq_c]; % deviation from equilibrium

% --- A1: Step ---  Dv = +0.2 around v=4.0
t_step_A  = 20;
dv_step_A = 0.2;
u_stepA   = @(t) v_c + dv_step_A*(t >= t_step_A);
du_stepA  = dv_step_A * (t_eval >= t_step_A);

[t_nl_sA, x_nl_sA] = ode45(...
    @(t,x) seir_ode(t,x,u_stepA(t),b_c,mu,beta_c,c_par,sigma,delta,epsilon,gamma,d_inf),...
    [0,t_end], x0_A, ode_opts);
[~,~,x_lin_sA] = lsim(sys_c, du_stepA, t_eval, dx0_A);

I_nl_sA  = interp1(t_nl_sA, x_nl_sA(:,3), t_eval, 'pchip');
I_lin_sA = x_lin_sA(:,3) + 0;   % I_eq = 0

rmse_sA  = sqrt(mean((I_nl_sA(mask) - I_lin_sA(mask)).^2));
nrmse_sA = rmse_sA / (max(I_nl_sA) - min(I_nl_sA) + eps);
fprintf('  Step (Dv=+%.1f around v=%.1f): RMSE=%.2e  NRMSE=%.4f\n',...
        dv_step_A, v_c, rmse_sA, nrmse_sA);

% --- A2: Sinusoidal --- Dv = 0.15*sin(0.1*t) around v=4.0
A_sinA   = 0.15;
om_sinA  = 0.1;
u_sinA   = @(t) v_c + A_sinA*sin(om_sinA*t);
du_sinA  = A_sinA * sin(om_sinA * t_eval);

[t_nl_sinA, x_nl_sinA] = ode45(...
    @(t,x) seir_ode(t,x,u_sinA(t),b_c,mu,beta_c,c_par,sigma,delta,epsilon,gamma,d_inf),...
    [0,t_end], x0_A, ode_opts);
[~,~,x_lin_sinA] = lsim(sys_c, du_sinA, t_eval, dx0_A);

I_nl_sinA  = interp1(t_nl_sinA, x_nl_sinA(:,3), t_eval, 'pchip');
I_lin_sinA = x_lin_sinA(:,3);

rmse_sinA  = sqrt(mean((I_nl_sinA(mask) - I_lin_sinA(mask)).^2));
nrmse_sinA = rmse_sinA / (max(I_nl_sinA) - min(I_nl_sinA) + eps);
fprintf('  Sin  (A=%.2f, w=%.1f around v=%.1f): RMSE=%.2e  NRMSE=%.4f\n\n',...
        A_sinA, om_sinA, v_c, rmse_sinA, nrmse_sinA);

% ============================================================
%  SCENARIO B: Groupmate's exact parameters (v=0.1)
%  Both nonlinear and linear use sys1 (built at v=0.1)
% ============================================================
fprintf('  --- Scenario B: Groupmate params (v=0.1, R0=5.6 > 1) ---\n');

% Correct initial conditions: groupmate's exact values
x0_B  = [0.95*S0_1; 2; 5; 0];
dx0_B = x0_B - [S0_1; 0; 0; 0];  % deviation from Set1 DFE

% --- B1: Step --- v: 0.1 -> 0.2 (deviation +0.1 around v=0.1)
v_before  = 0.1;
v_after   = 0.2;
t_step_B  = 20;
dv_step_B = v_after - v_before;
u_stepB   = @(t) v_before + dv_step_B*(t >= t_step_B);
du_stepB  = dv_step_B * (t_eval >= t_step_B);

[t_nl_sB, x_nl_sB] = ode45(...
    @(t,x) seir_ode(t,x,u_stepB(t),b1,mu,beta1,c_par,sigma,delta,epsilon,gamma,d_inf),...
    [0,t_end], x0_B, ode_opts);
[~,~,x_lin_sB] = lsim(sys1, du_stepB, t_eval, dx0_B);  % sys1 at v=0.1

I_nl_sB  = interp1(t_nl_sB, x_nl_sB(:,3), t_eval, 'pchip');
I_lin_sB = x_lin_sB(:,3) + 0;

rmse_sB  = sqrt(mean((I_nl_sB(mask) - I_lin_sB(mask)).^2));
nrmse_sB = rmse_sB / (max(I_nl_sB) - min(I_nl_sB) + eps);
fprintf('  Step (v:0.1->0.2): RMSE=%.2e  NRMSE=%.4f\n', rmse_sB, nrmse_sB);

% --- B2: Sinusoidal --- v(t)=0.1+0.05*sin(0.02*t) (groupmate's exact)
v0_sin  = 0.1;
A_sin   = 0.05;
om_sin  = 0.02;
u_sinB  = @(t) v0_sin + A_sin*sin(om_sin*t);
du_sinB = A_sin * sin(om_sin * t_eval);

[t_nl_sinB, x_nl_sinB] = ode45(...
    @(t,x) seir_ode(t,x,u_sinB(t),b1,mu,beta1,c_par,sigma,delta,epsilon,gamma,d_inf),...
    [0,t_end], x0_B, ode_opts);
[~,~,x_lin_sinB] = lsim(sys1, du_sinB, t_eval, dx0_B);  % sys1 at v=0.1

I_nl_sinB  = interp1(t_nl_sinB, x_nl_sinB(:,3), t_eval, 'pchip');
I_lin_sinB = x_lin_sinB(:,3);

rmse_sinB  = sqrt(mean((I_nl_sinB(mask) - I_lin_sinB(mask)).^2));
nrmse_sinB = rmse_sinB / (max(I_nl_sinB) - min(I_nl_sinB) + eps);
fprintf('  Sin  (v=0.1+0.05*sin(0.02t)): RMSE=%.2e  NRMSE=%.4f\n\n',...
        rmse_sinB, nrmse_sinB);
fprintf('  Validity range: Linearized model valid for approx +/-5%% from equilibrium\n\n');

% Alias for figures (use Scenario A for clean RMSE plots, B for groupmate comparison)
t_nl_s  = t_nl_sA;   x_nl_s  = x_nl_sA;
I_nl_s  = I_nl_sA;   I_lin_s = I_lin_sA;
rmse_s  = rmse_sA;   nrmse_s = nrmse_sA;
t_step  = t_step_A;  v_before_label = v_c; v_after_label = v_c + dv_step_A;

t_nl_sin  = t_nl_sinB;  % Fig 4/9 use Scenario B (groupmate's exact for visual)
I_nl_sin  = I_nl_sinB;
I_lin_sin = I_lin_sinB;
rmse_sin  = rmse_sinB;  nrmse_sin = nrmse_sinB;

%% =========================================================
%  SEC 12: g1 vs g2 COMPARISON SIMULATION
% =========================================================
fprintf('[SEC 12] g1 vs g2 comparison simulation...\n');

% Both use Scenario A's consistent setup (v=4.0)
u_g_compare = @(t) v_c + A_sinA*sin(om_sinA*t);

[t_g1_raw, x_g1_raw] = ode45(...
    @(t,x) seir_ode(t,x,u_g_compare(t),b_c,mu,beta_c,c_par,sigma,delta,epsilon,gamma,d_inf),...
    [0,t_end], x0_A, ode_opts);

% g2: set c_par=0 (bilinear -- no saturation denominator)
[t_g2, x_g2] = ode45(...
    @(t,x) seir_ode(t,x,u_g_compare(t),b_c,mu,beta_c,0,sigma,delta,epsilon,gamma,d_inf),...
    [0,t_end], x0_A, ode_opts);

I_g1 = interp1(t_g1_raw, x_g1_raw(:,3), t_eval, 'pchip');
I_g2_raw = interp1(t_g2, x_g2(:,3), t_eval, 'pchip');
rmse_g1g2 = sqrt(mean((I_g1 - I_g2_raw).^2));
fprintf('  RMSE between g1 and g2 in I(t) = %.6f\n', rmse_g1g2);
fprintf('  At I=0 (DFE), g1 and g2 give identical linearization\n');
fprintf('  At large I, g1 saturates while g2 grows unbounded\n\n');

%% =========================================================
%  SEC 13: ALL 9 FIGURES
% =========================================================
fprintf('[SEC 13] Generating 9 figures...\n\n');

% ---- Figure 1: Pole-Zero Map (custom zoomed) ----
fh1 = figure('Name','Fig 1 - Pole-Zero Map','NumberTitle','off',...
             'Position',[50,50,820,520]);
p_re = real(sys_poles_c);
p_im = imag(sys_poles_c);
z_re = real(sys_zeros_c);
z_im = imag(sys_zeros_c);
plot(p_re, p_im, 'rx', 'MarkerSize',14, 'LineWidth',2.5, 'DisplayName','Poles');
hold on;
if ~isempty(sys_zeros_c)
    plot(z_re, z_im, 'bo', 'MarkerSize',12, 'LineWidth',2.5, ...
         'DisplayName','Zeros', 'MarkerFaceColor','none');
end
xline(0,'k--','LineWidth',1.5,'DisplayName','Stability boundary (Re=0)');
for k = 1:length(sys_poles_c)
    text(p_re(k), p_im(k)+0.002, ...
         sprintf('  \\lambda_{%d} = %.4f', k, p_re(k)), ...
         'FontSize',10, 'Color',[0.8 0 0]);
end
if ~isempty(sys_zeros_c)
    text(z_re(1), z_im(1)+0.002, ...
         sprintf('  z = %.4f', z_re(1)), 'FontSize',10, 'Color',[0 0 0.8]);
end
all_re_pz = [p_re; z_re];
margin_pz = max(abs(all_re_pz)) * 0.2 + 0.05;
xlim([min(all_re_pz) - margin_pz*5, max(all_re_pz) + margin_pz*2]);
ylim([-0.06, 0.06]);
xlabel('Real Axis', 'FontSize',12);
ylabel('Imaginary Axis', 'FontSize',12);
title({'Pole-Zero Map of Linearized SEIR System  (y = I, Control case v_0=4.0)',...
       'All poles on real axis \Rightarrow purely exponential response, no oscillation'},...
      'FontSize',12,'FontWeight','bold');
legend('Location','best','FontSize',10); grid on;
annotation('textbox',[0.55 0.08 0.38 0.13],...
    'String',{sprintf('R_0 = %.4f < 1  \\Rightarrow  Stable', R0_c),...
               'All poles: Re(\lambda) < 0'},...
    'FontSize',11,'BackgroundColor',[0.9 1.0 0.9],...
    'EdgeColor',[0 0.6 0],'FontWeight','bold');
savefig_local(fh1, SAVE_FIGS, OUTPUT_DIR, 'fig1_pole_zero_map.png');

% ---- Figure 2: Sensitivity Analysis ----
fh2 = figure('Name','Fig 2 - Sensitivity','NumberTitle','off',...
             'Position',[80,50,920,480]);
bv_axis = beta_c * beta_mults;
subplot(1,2,1);
plot(bv_axis, all_eigs_sens', 'o-', 'LineWidth',2, 'MarkerSize',7);
hold on;
yline(0,'k--','LineWidth',1.8,'DisplayName','Boundary (Re=0)');
xlabel('\beta (Transmission Rate)','FontSize',12);
ylabel('Re(\lambda_i)','FontSize',12);
title('Eigenvalues vs beta  (+/-5%, +/-10%)','FontSize',12,'FontWeight','bold');
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','Boundary',...
       'Location','best','FontSize',9);
xticks(bv_axis);
xticklabels({'-10%','-5%','0%','+5%','+10%'});
grid on;
subplot(1,2,2);
plot(bv_axis, R0_sens, 'ms-', 'LineWidth',2, 'MarkerSize',8);
hold on;
yline(1,'r--','LineWidth',1.8,'DisplayName','R_0=1 threshold');
yline(R0_c,'g:','LineWidth',1.5,'DisplayName',sprintf('Nominal R_0=%.4f',R0_c));
xlabel('\beta (Transmission Rate)','FontSize',12);
ylabel('R_0','FontSize',12);
title('R_0 vs \beta Sensitivity','FontSize',12,'FontWeight','bold');
legend('Location','best','FontSize',9); grid on;
xticks(bv_axis);
xticklabels({'-10%','-5%','0%','+5%','+10%'});
sgtitle('Figure 2 — Sensitivity Analysis: beta +/-10%','FontSize',13,'FontWeight','bold');
savefig_local(fh2, SAVE_FIGS, OUTPUT_DIR, 'fig2_sensitivity.png');

% ---- Figure 3: Step Response ----
fh3 = figure('Name','Fig 3 - Step Response','NumberTitle','off',...
             'Position',[100,50,1000,620]);
subplot(2,1,1);
plot(t_eval, I_nl_s, 'b-',  'LineWidth',2, 'DisplayName','Nonlinear I(t)  [ode45]');
hold on;
plot(t_eval, I_lin_s,'r--', 'LineWidth',2, 'DisplayName','Linearized I(t)  [lsim]');
xline(t_step_A,'k:','LineWidth',1.5,'DisplayName',sprintf('Step at t=%d days',t_step_A));
ylabel('Infected  I(t)  [persons]','FontSize',12);
title(sprintf('Figure 3 — Step Input (Dv=+%.1f at t=%d, v0=%.1f): Nonlinear vs Linearized',...
              dv_step_A, t_step_A, v_c),'FontSize',12,'FontWeight','bold');
legend('Location','best','FontSize',10); xlim([0,t_end]); grid on;
subplot(2,1,2);
plot(t_eval, abs(I_nl_s - I_lin_s),'Color',[0.9 0.4 0.0],'LineWidth',1.8);
xlabel('Time (days)','FontSize',12);
ylabel('|Error|  [persons]','FontSize',12);
title(sprintf('Absolute Error  |  RMSE = %.2e  |  NRMSE = %.4f  (Scenario A, v=4.0)', rmse_s, nrmse_s),'FontSize',11);
xlim([0,t_end]); grid on;
savefig_local(fh3, SAVE_FIGS, OUTPUT_DIR, 'fig3_step_response.png');

% ---- Figure 4: Sinusoidal Response ----
fh4 = figure('Name','Fig 4 - Sinusoidal','NumberTitle','off',...
             'Position',[120,50,1000,620]);
subplot(2,1,1);
plot(t_eval, I_nl_sin, 'b-',  'LineWidth',2, 'DisplayName','Nonlinear I(t)  [ode45]');
hold on;
plot(t_eval, I_lin_sin,'r--', 'LineWidth',2, 'DisplayName','Linearized I(t)  [lsim]');
ylabel('Infected  I(t)  [persons]','FontSize',12);
title('Figure 4 — Sinusoidal v(t)=0.1+0.05sin(0.02t): Nonlinear vs Linearized (Scenario B, v=0.1)','FontSize',12,'FontWeight','bold');
legend('Location','best','FontSize',10); xlim([0,t_end]); grid on;
subplot(2,1,2);
plot(t_eval, abs(I_nl_sin - I_lin_sin),'Color',[0.55 0 0.75],'LineWidth',1.8);
xlabel('Time (days)','FontSize',12);
ylabel('|Error|  [persons]','FontSize',12);
title(sprintf('Absolute Error  |  RMSE=%.2e  NRMSE=%.4f  (Scenario B)', rmse_sin, nrmse_sin),'FontSize',11);
xlim([0,t_end]); grid on;
savefig_local(fh4, SAVE_FIGS, OUTPUT_DIR, 'fig4_sinusoidal.png');

% ---- Figure 5: Full SEIR State Trajectories ----
fh5 = figure('Name','Fig 5 - SEIR States','NumberTitle','off',...
             'Position',[140,50,1100,720]);
state_names = {'S(t)  Susceptible','E(t)  Exposed','I(t)  Infected','R(t)  Recovered'};
R0_eq_c_fig = v_c * S0_c / mu;
eq_vals_c   = [S0_c, 0, 0, R0_eq_c_fig];
clrs_plot   = {'b','r','m','g'};
for k = 1:4
    subplot(2,2,k);
    plot(t_nl_s, x_nl_s(:,k), '-','Color',clrs_plot{k},'LineWidth',2,'DisplayName','Trajectory');
    hold on;
    yline(eq_vals_c(k),'k--','LineWidth',1.3,'DisplayName','Equilibrium');
    xlabel('Time (days)','FontSize',11);
    ylabel('Population','FontSize',11);
    title(state_names{k},'FontSize',12,'FontWeight','bold');
    legend('Location','best','FontSize',9); grid on;
end
sgtitle('Figure 5 — SEIR Nonlinear Trajectories  (Step v: 0.1 to 0.2, Control case R0 < 1)',...
        'FontSize',12,'FontWeight','bold');
savefig_local(fh5, SAVE_FIGS, OUTPUT_DIR, 'fig5_seir_states.png');

% ---- Figure 6: Structural Analysis Singular Values ----
fh6 = figure('Name','Fig 6 - Structural','NumberTitle','off',...
             'Position',[160,50,900,420]);
subplot(1,2,1);
bar(svd(Co_mat),'FaceColor','b','EdgeColor','k');
title(sprintf('Controllability Matrix \\sigma_i  (rank=%d/%d)',rank_Co,n),...
      'FontSize',12,'FontWeight','bold');
xlabel('Index','FontSize',11); ylabel('\sigma_i','FontSize',11); grid on;
subplot(1,2,2);
bar(svd(Ob_mat),'FaceColor','r','EdgeColor','k');
title(sprintf('Observability Matrix \\sigma_i  (rank=%d/%d)',rank_Ob,n),...
      'FontSize',12,'FontWeight','bold');
xlabel('Index','FontSize',11); ylabel('\sigma_i','FontSize',11); grid on;
sgtitle('Figure 6 — SVD Structural Analysis: Singular Values of Co and Ob',...
        'FontSize',12,'FontWeight','bold');
savefig_local(fh6, SAVE_FIGS, OUTPUT_DIR, 'fig6_structural.png');

% ---- Figure 7: Two Parameter Sets Comparison ----
fh7 = figure('Name','Fig 7 - Two Sets','NumberTitle','off',...
             'Position',[180,50,1100,720]);
x0_1 = [0.95*S0_1; 2; 5; 0];
x0_2 = [0.95*S0_2; 2; 5; 0];
t2   = linspace(0, 200, 3000);
[t_s1, x_s1] = ode45(@(t,x) seir_ode(t,x,v1,b1,mu,beta1,c_par,sigma,delta,epsilon,gamma,d_inf),...
    [0,200], x0_1, ode_opts);
[t_s2, x_s2] = ode45(@(t,x) seir_ode(t,x,v2,b2,mu,beta2,c_par,sigma,delta,epsilon,gamma,d_inf),...
    [0,200], x0_2, ode_opts);
names7 = {'S(t)','E(t)','I(t)','R(t)'};
for k = 1:4
    subplot(2,2,k);
    plot(t_s1, x_s1(:,k),'b-','LineWidth',2,...
         'DisplayName',sprintf('Set 1 (v=%.1f, R_0=%.2f)', v1, R0_1));
    hold on;
    try
        plot(t_s2, x_s2(:,k),'r--','LineWidth',2,...
             'DisplayName',sprintf('Set 2 (v=%.0f, R_0=%.0f)', v2, R0_2));
    catch
    end
    xlabel('Time (days)','FontSize',11);
    ylabel('Population','FontSize',11);
    title(names7{k},'FontSize',12,'FontWeight','bold');
    legend('Location','best','FontSize',8); grid on;
end
sgtitle('Figure 7 — Two Parameter Sets: Set 1 (v=0.1, beta=0.01) vs Set 2 (v=0, beta=0.9)',...
        'FontSize',12,'FontWeight','bold');
savefig_local(fh7, SAVE_FIGS, OUTPUT_DIR, 'fig7_two_sets.png');

% ---- Figure 8: g1 vs g2 Comparison ----
fh8 = figure('Name','Fig 8 - g1 vs g2','NumberTitle','off',...
             'Position',[200,50,1000,480]);
subplot(1,2,1);
plot(t_eval, I_g1,'b-','LineWidth',2,'DisplayName','g1 = betaSI/(1+cI^2)');
hold on;
I_g2_plot = I_g2_raw;
plot(t_eval, I_g2_plot,'r--','LineWidth',2,'DisplayName','g2 = betaSI');
xlabel('Time (days)','FontSize',12); ylabel('Infected I(t)','FontSize',12);
title('Infected I(t): g_1 vs g_2','FontSize',12,'FontWeight','bold');
legend('Location','best','FontSize',10); grid on;
subplot(1,2,2);
diff_g1g2 = abs(I_g1 - I_g2_plot);
plot(t_eval, diff_g1g2,'g-','LineWidth',2);
fill([t_eval; flipud(t_eval)], [diff_g1g2; zeros(size(diff_g1g2))], ...
     'g', 'FaceAlpha',0.2, 'EdgeColor','none');
xlabel('Time (days)','FontSize',12);
ylabel('|g_1 - g_2| in I(t)','FontSize',12);
title(sprintf('|Difference|  (RMSE=%.4f)', rmse_g1g2),'FontSize',12,'FontWeight','bold');
grid on;
sgtitle('Figure 8 — Saturated g_1 vs Bilinear g_2 Incidence Functions',...
        'FontSize',12,'FontWeight','bold');
savefig_local(fh8, SAVE_FIGS, OUTPUT_DIR, 'fig8_g1_vs_g2.png');

% ---- Figure 9: Sinusoidal -- groupmate's exact parameters ----
fh9 = figure('Name','Fig 9 - Sinusoidal (groupmate params)','NumberTitle','off',...
             'Position',[220,50,1000,620]);
subplot(2,1,1);
plot(t_eval, I_nl_sin, 'b-',  'LineWidth',2, 'DisplayName','Nonlinear I(t)');
hold on;
plot(t_eval, I_lin_sin,'r--', 'LineWidth',2, 'DisplayName','Linearized I(t)');
ylabel('Infected  I(t)  [persons]','FontSize',12);
title(sprintf('Figure 9 — Sinusoidal v(t)=0.1+0.05sin(0.02t): RMSE=%.2e, NRMSE=%.4f (Scenario B, sys1 at v=0.1)',...
              rmse_sinB, nrmse_sinB),...
      'FontSize',11,'FontWeight','bold');
legend('Location','best','FontSize',10); xlim([0,t_end]); grid on;
subplot(2,1,2);
plot(t_eval, abs(I_nl_sin - I_lin_sin),'Color',[0.55 0 0.75],'LineWidth',1.8);
xlabel('Time (days)','FontSize',12);
ylabel('|Error|','FontSize',12);
title('Absolute Error — Validity range: approx ±5% from equilibrium','FontSize',11);
xlim([0,t_end]); grid on;
savefig_local(fh9, SAVE_FIGS, OUTPUT_DIR, 'fig9_sinusoidal_groupmate.png');

%% =========================================================
%  SEC 14: COMPLETE SUMMARY
% =========================================================
% Helper strings for summary
if R0_1 <= 1, s1 = 'Stable DFE'; else, s1 = 'UNSTABLE DFE'; end
if R0_2 <= 1, s2 = 'Stable DFE'; else, s2 = 'UNSTABLE DFE'; end
if R0_c <= 1, sc = 'STABLE DFE'; else, sc = 'UNSTABLE DFE'; end
if rank_Co == n, ctrl_str = 'Fully Controllable'; else, ctrl_str = 'PARTIALLY Controllable'; end
if rank_Ob == n, obs_str = 'Fully Observable'; else, obs_str = 'PARTIALLY Observable'; end

fprintf('=================================================================');
fprintf('                  PHASE 1 COMPLETE SUMMARY\n');
fprintf('=================================================================\n');
fprintf('  EQUILIBRIA:\n');
fprintf('    Set 1 (v=%.1f)  : S0=%.2f   R0=%.4f  --> %s\n', ...
        v1, S0_1, R0_1, s1);
fprintf('    Set 2 (v=%.0f)    : S0=%.2f  R0=%.2f --> %s\n', ...
        v2, S0_2, R0_2, s2);
fprintf('    Control (v=%.1f): S0=%.4f  R0=%.4f --> %s\n', ...
        v_c, S0_c, R0_c, sc);
fprintf('  ENDEMIC E** (Set 2): S**=%.2f, E**=%.2f, I**=%.2f, R**=%.2f\n',...
        S_star, E_star, I_star, R_star);
fprintf('  EIGENVALUES (Control v=%.1f):\n', v_c);
fprintf('    [%.5f, %.5f, %.5f, %.5f]\n',...
        eigs_c(1),eigs_c(2),eigs_c(3),eigs_c(4));
if all(eigs_c<0), stab_all='YES'; else, stab_all='NO'; end
fprintf('    All stable? %s\n', stab_all);
fprintf('  CANONICAL FORMS:\n');
fprintf('    Jordan : J = diag(%.4f, %.4f, %.4f, %.4f)\n',...
        J_jord(1,1),J_jord(2,2),J_jord(3,3),J_jord(4,4));
fprintf('    CCF    : Bc = [0;0;0;1]  last row = [-a0,-a1,-a2,-a3]\n');
fprintf('    OCF    : Co = [0 0 0 1]  last col = [-a0;-a1;-a2;-a3]\n');
fprintf('  STRUCTURAL ANALYSIS:\n');
fprintf('    rank(Co)  = %d / %d  --> %s (all 5 methods)\n',...
        rank_Co, n, ctrl_str);
fprintf('    rank(Ob)  = %d / %d  --> %s (y=I only)\n',...
        rank_Ob, n, obs_str);
fprintf('    rank(Ob)  = %d / %d  --> Fully Observable (y=[S,E,I,R])\n', rank(Ob_full), n);
fprintf('    Minimal order: %d\n', n_min);
fprintf('  SIMULATION:\n');
fprintf('    Step     : RMSE=%.2e, NRMSE=%.4f\n', rmse_s, nrmse_s);
fprintf('    Sinusoidal: RMSE=%.2e, NRMSE=%.4f\n', rmse_sin, nrmse_sin);
fprintf('    g1 vs g2 difference RMSE: %.6f\n', rmse_g1g2);
fprintf('    Valid range: ~+/-5%% from equilibrium\n');
fprintf('  SENSITIVITY: beta +/-10%% --> eigenvalues nearly unchanged\n');
fprintf('=================================================================\n\n');

if SAVE_FIGS
    fprintf('All 9 figures saved to: %s\n', OUTPUT_DIR);
else
    fprintf('All 9 figures displayed on screen.\n');
    fprintf('Set SAVE_FIGS = true to save as PNG files.\n');
end
fprintf('DONE.\n');


%% =========================================================
%  LOCAL FUNCTIONS  (must be at bottom of script file)
% =========================================================

% ---------------------------------------------------------
% seir_ode: Right-hand side of SEIR system (works for g1 and g2)
% Set c_par=0 for bilinear (g2) behaviour
% ---------------------------------------------------------
function dxdt = seir_ode(~, x, u, b, mu, beta, c_par, ...
                          sigma, delta, epsilon, gamma, d_inf)
    S = x(1);  E = x(2);  I = x(3);  R = x(4);
    g    = beta * S * I / (1 + c_par * I^2);   % g1 (set c_par=0 for g2)
    dxdt = [ b  - g - (mu + u)*S;
             g  - (epsilon + sigma + mu + delta)*E;
             sigma*E - (mu + gamma + d_inf)*I;
             u*S + gamma*I + epsilon*E - mu*R ];
end

% ---------------------------------------------------------
% linearize_SEIR: Numerical central-difference Jacobian
% Works for both g1 and g2 (controlled via c_par)
% ---------------------------------------------------------
function [A, B, C, D] = linearize_SEIR(x_eq, beta, v, mu, ...
                          epsilon, sigma, delta, gamma, d_inf)
    S0 = x_eq(1);
    % At DFE (I=0): dg/dS = 0, dg/dI = beta*S0  (same for g1 and g2)
    dg_dS = 0;
    dg_dI = beta * S0;

    A = [ -(dg_dS) - (mu+v),     0,                         -(dg_dI),          0;
           dg_dS,              -(epsilon+sigma+mu+delta),     dg_dI,             0;
           0,                   sigma,                      -(mu+gamma+d_inf),  0;
           v,                   epsilon,                     gamma,            -mu];

    B = [-S0; 0; 0; S0];
    C = [0, 0, 1, 0];
    D = 0;
end

% ---------------------------------------------------------
% compute_K: Used by fzero to find I** (endemic equilibrium)
% ---------------------------------------------------------
function val = compute_K(I, b, mu, v, beta, c, sigma, delta, epsilon, gamma, d)
    S_num = b/(mu+v) - (sigma+mu+delta+epsilon)*(mu+gamma+d) / (sigma*(mu+v)) * I;
    if S_num <= 0
        val = -1;
        return;
    end
    g_val = beta * S_num * I / (1 + c * I^2);
    val = g_val / I - (sigma+mu+delta+epsilon)*(mu+gamma+d) / sigma;
end

% ---------------------------------------------------------
% savefig_local: Save figure if SAVE_FIGS is true
% ---------------------------------------------------------
function savefig_local(fig, save_flag, out_dir, fname)
    if save_flag
        saveas(fig, fullfile(out_dir, fname));
        fprintf('  Saved: %s\n', fname);
    end
end
