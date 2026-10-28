```matlab
%% SCRI 2026/27 - TP1: Estimacao de Estado (WLS)
%  Modelo, Ybus, vetor de estado, vetor de medidas e funcoes de medida
%  h(x) escritas em funcao de Y = G + jB (calculo de h na iteracao 0)
clear; clc;


% ==================== MODELO ====================
% --- Admitancias serie das linhas (pu) ---
y12a = 2 - 10j;          % linha a entre 1 e 2
y12b = 2 - 10j;          % linha b entre 1 e 2 (paralela)
y13  = 2 - 10j;
y23  = 2 - 10j;


% --- Admitancia shunt do banco de condensadores (barramento 1) ---
Bcap = 2.2;              % shunt = j*Bcap  ->  Qgbat = Bcap*V1^2


% --- Medicoes (pu) ---
V1    = 1.10;
V2    = 1.00;
V3    = 1.07;
P21a  = -0.79;
P21b  =  0.79;
P31   =  1.71;
P32   =  2.34;
Q12a  =  0.97;
Q32   =  1.22;
Qgbat =  2.66;
P1v   =  0;              % medida virtual: barramento 1 sem carga/geracao
Q1v   =  0;              % medida virtual


% --- Variancias dos erros de medida ---
varV1    = 0.0002;
varV2    = 0.0002;
varV3    = 0.0002;
varP21a  = 0.001;
varP21b  = 0.001;
varP31   = 0.001;
varP32   = 0.001;
varQ12a  = 0.005;
varQ32   = 0.005;
varQgbat = 0.005;
varP1v   = 0.00001;
varQ1v   = 0.00001;


% --- Angulos ja conhecidos ---
t2 = 0;                  % barramento 2 = slack (referencia). t1 e t3 sao estimados.


% --- Parametros do metodo ---
tol   = 1e-4;            % tolerancia em max|dx|
maxit = 50;
alfa  = 0.05;            % nivel de significancia do teste chi^2
limiar_rN = 3;           % limiar para |rN|


% Estrutura com os parametros da rede (passada as funcoes)
rede.g12a = real(y12a);  rede.b12a = imag(y12a);
rede.g12b = real(y12b);  rede.b12b = imag(y12b);
rede.g13  = real(y13);   rede.b13  = imag(y13);
rede.g23  = real(y23);   rede.b23  = imag(y23);
rede.Bcap = Bcap;        rede.t2   = t2;


% ================== MATRIZ Ybus ==================
Y11 =  y12a + y12b + y13 + 1j*Bcap;
Y22 =  y12a + y12b + y23;
Y33 =  y13 + y23;
Y12 = -(y12a + y12b);    Y21 = Y12;
Y13 = -y13;              Y31 = Y13;
Y23 = -y23;              Y32 = Y23;
Ybus = [Y11 Y12 Y13;
        Y21 Y22 Y23;
        Y31 Y32 Y33];


BusNames = {'ONE', 'TWO', 'THREE'};
disp_Ybus = array2table(Ybus, 'RowNames', BusNames, 'VariableNames', BusNames)
```

| |ONE|TWO|THREE|
|:--:|:--:|:--:|:--:|
|1 ONE|6.0000 -27.8000i|-4.0000 +20.0000i|-2.0000 +10.0000i|
|2 TWO|-4.0000 +20.0000i|6.0000 -30.0000i|-2.0000 +10.0000i|
|3 THREE|-2.0000 +10.0000i|-2.0000 +10.0000i|4.0000 -20.0000i|

```matlab


% ================ VETOR DE ESTADO ================
% x = [t1, t3, V1, V2, V3]'   (angulos em rad)
StateNames  = {'t1', 't3', 'V1', 'V2', 'V3'};
StateNamesD = {'t1_graus', 't3_graus', 'V1_pu', 'V2_pu', 'V3_pu'};   % para as tabelas
n = numel(StateNames);


% ================ VETOR DE MEDIDAS ================
MeasNames = {'V1','V2','V3','P21a','P21b','P31','P32','Q12a','Q32','Qgbat','P1v','Q1v'};
z   = [V1; V2; V3; P21a; P21b; P31; P32; Q12a; Q32; Qgbat; P1v; Q1v];
var = [varV1; varV2; varV3; varP21a; varP21b; varP31; varP32; ...
       varQ12a; varQ32; varQgbat; varP1v; varQ1v];
virtual = ismember(MeasNames, {'P1v','Q1v'})';


MeasCols = {'z', 'variancia'};
disp_medidas = array2table([z, var], 'RowNames', MeasNames, 'VariableNames', MeasCols)
```

| |z|variancia|
|:--:|:--:|:--:|
|1 V1|1.1000|2.0000e-04|
|2 V2|1|2.0000e-04|
|3 V3|1.0700|2.0000e-04|
|4 P21a|-0.7900|1.0000e-03|
|5 P21b|0.7900|1.0000e-03|
|6 P31|1.7100|1.0000e-03|
|7 P32|2.3400|1.0000e-03|
|8 Q12a|0.9700|0.0050|
|9 Q32|1.2200|0.0050|
|10 Qgbat|2.6600|0.0050|
|11 P1v|0|1.0000e-05|
|12 Q1v|0|1.0000e-05|

```matlab


% ================= MATRIZ DE PESOS W =================
% W = R^-1, com R = diag(variancias dos erros de medida)
R = diag(var);
W = diag(1 ./ var);
disp_W = array2table(diag(W), 'RowNames', MeasNames, 'VariableNames', {'peso_1_sobre_var'})
```

| |peso_1_sobre_var|
|:--:|:--:|
|1 V1|5000|
|2 V2|5000|
|3 V3|5000|
|4 P21a|1000|
|5 P21b|1000|
|6 P31|1000|
|7 P32|1000|
|8 Q12a|200|
|9 Q32|200|
|10 Qgbat|200|
|11 P1v|1.0000e+05|
|12 Q1v|1.0000e+05|

```matlab


% ==================================================================
%  FUNCOES DE MEDIDA h(x) E JACOBIANA H(x) EM FUNCAO DE Y = G + jB
% ==================================================================
% Convencao: G_ik + jB_ik = Y_ik = -y_ik (elemento fora da diagonal da Ybus)
%   P_ik = Vi*Vk*(G_ik*cos(t_ik) + B_ik*sin(t_ik)) - G_ik*Vi^2
%   Q_ik = Vi*Vk*(G_ik*sin(t_ik) - B_ik*cos(t_ik)) + B_ik*Vi^2   (b_shunt,linha = 0)
% Jacobiana: colunas 1 = d/dt1 | 2 = d/dt3 | 3 = d/dV1 | 4 = d/dV2 | 5 = d/dV3
% Derivadas de um fluxo de i para k, com A_ik = G*cos(t_ik) + B*sin(t_ik)
%                                       e D_ik = G*sin(t_ik) - B*cos(t_ik):
%   dP/dti = -Vi*Vk*D_ik     dP/dtk = +Vi*Vk*D_ik
%   dP/dVi =  Vk*A_ik - 2*G*Vi        dP/dVk =  Vi*A_ik
%   dQ/dti = +Vi*Vk*A_ik     dQ/dtk = -Vi*Vk*A_ik
%   dQ/dVi =  Vk*D_ik + 2*B*Vi        dQ/dVk =  Vi*D_ik


% --- G e B de cada LINHA individual (para os fluxos) ---
G12a = -real(y12a);  B12a = -imag(y12a);
G12b = -real(y12b);  B12b = -imag(y12b);
G13  = -real(y13);   B13  = -imag(y13);
G23  = -real(y23);   B23  = -imag(y23);


% --- G e B da Ybus COMPLETA (para as medidas virtuais e Qgbat) ---
GY = real(Ybus);     BY = imag(Ybus);


% ==================================================================
%  PROCESSO ITERATIVO (GAUSS-NEWTON) - ATE max|dx| < tol
% ==================================================================
m = numel(z);                      % numero de medidas
k = m - n;                         % graus de liberdade do teste chi^2
if exist('chi2inv', 'file')
    chi2_crit = chi2inv(1 - alfa, k);
else                               % tabela chi^2 a 95% (k = 1..10)
    tabela_chi2 = [3.84 5.99 7.81 9.49 11.07 12.59 14.07 15.51 16.92 18.31];
    chi2_crit = tabela_chi2(k);
end


x  = [0; 0; 1; 1; 1];              % flat start: x0 = [t1 t3 V1 V2 V3]'
dx = Inf(n, 1);                    % correcao anterior (Inf para nao parar na 1a avaliacao)
it = 0;                            % contador de iteracoes
X_hist = x;  dx_hist = [];         % historico (opcional)


while true
    % ---------- 1) zest = h(x) e H(x) no estado atual x ----------
    t1 = x(1);  t3 = x(2);  V1 = x(3);  V2 = x(4);  V3 = x(5);


    % --- Diferencas angulares ---
    t12 = t1 - t2;   t21 = t2 - t1;
    t31 = t3 - t1;   t13 = t1 - t3;
    t32 = t3 - t2;


    % --- Funcoes de medida ---
    h1  = V1;                                                                  % V1
    h2  = V2;                                                                  % V2
    h3  = V3;                                                                  % V3
    h4  = V2*V1*(G12a*cos(t21) + B12a*sin(t21)) - G12a*V2^2;                   % P21a
    h5  = V2*V1*(G12b*cos(t21) + B12b*sin(t21)) - G12b*V2^2;                   % P21b
    h6  = V3*V1*(G13 *cos(t31) + B13 *sin(t31)) - G13 *V3^2;                   % P31
    h7  = V3*V2*(G23 *cos(t32) + B23 *sin(t32)) - G23 *V3^2;                   % P32
    h8  = V1*V2*(G12a*sin(t12) - B12a*cos(t12)) + B12a*V1^2;                   % Q12a
    h9  = V3*V2*(G23 *sin(t32) - B23 *cos(t32)) + B23 *V3^2;                   % Q32
    h10 = (BY(1,1) + BY(1,2) + BY(1,3))*V1^2;                                  % Qgbat = Bcap*V1^2
    h11 = V1*( V1*GY(1,1) ...                                                  % P1 (virtual)
             + V2*(GY(1,2)*cos(t12) + BY(1,2)*sin(t12)) ...
             + V3*(GY(1,3)*cos(t13) + BY(1,3)*sin(t13)) );
    h12 = V1*( -V1*BY(1,1) ...                                                 % Q1 (virtual)
             + V2*(GY(1,2)*sin(t12) - BY(1,2)*cos(t12)) ...
             + V3*(GY(1,3)*sin(t13) - BY(1,3)*cos(t13)) );


    zest = [h1; h2; h3; h4; h5; h6; h7; h8; h9; h10; h11; h12];   % zest = h(x)


    % ---- Jacobiana H(x) ----
    % --- Termos A e D de cada linha individual (G, B da linha) ---
    A21  = G12a*cos(t21) + B12a*sin(t21);   D21  = G12a*sin(t21) - B12a*cos(t21);  % linhas a e b (identicas)
    A31  = G13 *cos(t31) + B13 *sin(t31);   D31  = G13 *sin(t31) - B13 *cos(t31);
    A32  = G23 *cos(t32) + B23 *sin(t32);   D32  = G23 *sin(t32) - B23 *cos(t32);
    A12a = G12a*cos(t12) + B12a*sin(t12);   D12a = G12a*sin(t12) - B12a*cos(t12);


    % --- Termos A e D da Ybus completa (para as medidas virtuais) ---
    A12Y = GY(1,2)*cos(t12) + BY(1,2)*sin(t12);   D12Y = GY(1,2)*sin(t12) - BY(1,2)*cos(t12);
    A13Y = GY(1,3)*cos(t13) + BY(1,3)*sin(t13);   D13Y = GY(1,3)*sin(t13) - BY(1,3)*cos(t13);


    H = zeros(12, 5);


    % --- V1, V2, V3 ---
    H(1,3) = 1;
    H(2,4) = 1;
    H(3,5) = 1;


    % --- P21a  (i = 2, k = 1; t2 fixo) ---
    H(4,1) =  V2*V1*D21;
    H(4,2) =  0;
    H(4,3) =  V2*A21;
    H(4,4) =  V1*A21 - 2*G12a*V2;
    H(4,5) =  0;


    % --- P21b  (identica a P21a, com G12b, B12b) ---
    H(5,1) =  V2*V1*D21;
    H(5,2) =  0;
    H(5,3) =  V2*A21;
    H(5,4) =  V1*A21 - 2*G12b*V2;
    H(5,5) =  0;


    % --- P31  (i = 3, k = 1) ---
    H(6,1) =  V3*V1*D31;
    H(6,2) = -V3*V1*D31;
    H(6,3) =  V3*A31;
    H(6,4) =  0;
    H(6,5) =  V1*A31 - 2*G13*V3;


    % --- P32  (i = 3, k = 2; t2 fixo) ---
    H(7,1) =  0;
    H(7,2) = -V3*V2*D32;
    H(7,3) =  0;
    H(7,4) =  V3*A32;
    H(7,5) =  V2*A32 - 2*G23*V3;


    % --- Q12a  (i = 1, k = 2; t2 fixo) ---
    H(8,1) =  V1*V2*A12a;
    H(8,2) =  0;
    H(8,3) =  V2*D12a + 2*B12a*V1;
    H(8,4) =  V1*D12a;
    H(8,5) =  0;


    % --- Q32  (i = 3, k = 2; t2 fixo) ---
    H(9,1) =  0;
    H(9,2) =  V3*V2*A32;
    H(9,3) =  0;
    H(9,4) =  V3*D32;
    H(9,5) =  V2*D32 + 2*B23*V3;


    % --- Qgbat = (BY11 + BY12 + BY13)*V1^2 ---
    H(10,3) = 2*(BY(1,1) + BY(1,2) + BY(1,3))*V1;


    % --- P1 virtual: P1 = V1*[V1*GY11 + V2*A12Y + V3*A13Y] ---
    H(11,1) = -V1*V2*D12Y - V1*V3*D13Y;
    H(11,2) =  V1*V3*D13Y;
    H(11,3) =  2*V1*GY(1,1) + V2*A12Y + V3*A13Y;
    H(11,4) =  V1*A12Y;
    H(11,5) =  V1*A13Y;


    % --- Q1 virtual: Q1 = V1*[-V1*BY11 + V2*D12Y + V3*D13Y] ---
    H(12,1) =  V1*V2*A12Y + V1*V3*A13Y;
    H(12,2) = -V1*V3*A13Y;
    H(12,3) = -2*V1*BY(1,1) + V2*D12Y + V3*D13Y;
    H(12,4) =  V1*D12Y;
    H(12,5) =  V1*D13Y;


    % ---------- 2) residuos, ganho, J e residuos normalizados ----------
    r     = z - zest;                          % residuos
    G     = H' * W * H;                        % matriz de ganho
    J     = sum(r.^2 ./ var);                  % funcao objetivo J(x)
    Omega = R - H * (G \ H');                  % covariancia dos residuos
    rN    = r ./ sqrt(abs(diag(Omega)));       % residuos normalizados


    % ---------- 3) criterio de paragem: max|dx| < tol ----------
    if max(abs(dx)) < tol
        break
    end
    if it >= maxit
        warning('Numero maximo de iteracoes (%d) atingido sem convergencia.', maxit);
        break
    end


    % ---------- 4) correcao: dx = [H'WH]^(-1) H'W [z - zest] ----------
    dx = G \ (H' * W * r);                     % equivale a inv(G) * H' * W * (z - zest)
    x  = x + dx;                               % x(it+1) = x(it) + dx
    it = it + 1;
    X_hist  = [X_hist, x];   dx_hist = [dx_hist, max(abs(dx))];
    fprintf('Iteracao %d:  max|dx| = %.4e\n', it, max(abs(dx)));
end
```

```matlabTextOutput
Iteracao 1:  max|dx| = 1.8036e-01
Iteracao 2:  max|dx| = 1.6753e-02
Iteracao 3:  max|dx| = 2.3929e-04
Iteracao 4:  max|dx| = 1.6623e-05
```

```matlab


% ==================================================================
%  RESULTADOS DA ULTIMA ITERACAO (max|dx| < tol)
% ==================================================================
fprintf('\nConvergiu: max|dx| = %.4e < tol = %.1e  (%d iteracoes)\n\n', max(abs(dx)), tol, it);
```

```matlabTextOutput
Convergiu: max|dx| = 1.6623e-05 < tol = 1.0e-04  (4 iteracoes)
```

```matlab


idx_reais = find(~virtual);                    % medidas reais (exclui P1v, Q1v)
[rN_max, q] = max(abs(rN(idx_reais)));
nome_max = MeasNames{idx_reais(q)};


Resumo_nomes = {'theta1 (graus)'; 'theta3 (graus)'; 'V1 (pu)'; 'V2 (pu)'; 'V3 (pu)'; ...
                'm / k = m - n'; 'Iteracoes'; 'J(x)'; 'chi2_{k;0,05}'; 'max |rN| (medidas reais)'};
Resumo_vals  = {rad2deg(x(1)); rad2deg(x(2)); x(3); x(4); x(5); ...
                sprintf('%d / %d', m, k); it; J; chi2_crit; ...
                sprintf('%.2f (%s)', rN_max, nome_max)};
disp_resumo = table(Resumo_vals, 'RowNames', Resumo_nomes, 'VariableNames', {'Valor'})
```

| |Valor|
|:--:|:--:|
|1 theta1 (graus)|2.5042|
|2 theta3 (graus)|9.3786|
|3 V1 (pu)|1.0890|
|4 V2 (pu)|0.9735|
|5 V3 (pu)|1.0995|
|6 m / k = m - n|'12 / 7'|
|7 Iteracoes|4|
|8 J(x)|2.3501e+03|
|9 chi2_{k;0,05}|14.0671|
|10 max \|rN\| (medidas reais)|'48.32 (P21b)'|

```matlab


% ==================================================================
%  2.a ESTIMACAO: DADOS INICIAIS, REMOVENDO P21b
% ==================================================================
ativas2 = true(numel(z), 1);
ativas2(strcmp(MeasNames, 'P21b')) = false;    % medida removida: P21b


Names2   = MeasNames(ativas2);
z_a      = z(ativas2);
var_a    = var(ativas2);
virtual2 = virtual(ativas2);
W2 = diag(1 ./ var_a);                         % pesos das medidas mantidas
R2 = diag(var_a);
m2 = nnz(ativas2);                             % numero de medidas utilizadas
k2 = m2 - n;                                   % graus de liberdade
if exist('chi2inv', 'file')
    chi2_crit2 = chi2inv(1 - alfa, k2);
else
    chi2_crit2 = tabela_chi2(k2);
end


x2  = [0; 0; 1; 1; 1];                         % flat start
dx2 = Inf(n, 1);
it2 = 0;


while true
    % ---------- 1) zest = h(x) e H(x) no estado atual x2 (todas as 12 linhas) ----------
    t1 = x2(1);  t3 = x2(2);  V1 = x2(3);  V2 = x2(4);  V3 = x2(5);


    % --- Diferencas angulares ---
    t12 = t1 - t2;   t21 = t2 - t1;
    t31 = t3 - t1;   t13 = t1 - t3;
    t32 = t3 - t2;


    % --- Funcoes de medida ---
    h1  = V1;                                                                  % V1
    h2  = V2;                                                                  % V2
    h3  = V3;                                                                  % V3
    h4  = V2*V1*(G12a*cos(t21) + B12a*sin(t21)) - G12a*V2^2;                   % P21a
    h5  = V2*V1*(G12b*cos(t21) + B12b*sin(t21)) - G12b*V2^2;                   % P21b
    h6  = V3*V1*(G13 *cos(t31) + B13 *sin(t31)) - G13 *V3^2;                   % P31
    h7  = V3*V2*(G23 *cos(t32) + B23 *sin(t32)) - G23 *V3^2;                   % P32
    h8  = V1*V2*(G12a*sin(t12) - B12a*cos(t12)) + B12a*V1^2;                   % Q12a
    h9  = V3*V2*(G23 *sin(t32) - B23 *cos(t32)) + B23 *V3^2;                   % Q32
    h10 = (BY(1,1) + BY(1,2) + BY(1,3))*V1^2;                                  % Qgbat = Bcap*V1^2
    h11 = V1*( V1*GY(1,1) ...                                                  % P1 (virtual)
             + V2*(GY(1,2)*cos(t12) + BY(1,2)*sin(t12)) ...
             + V3*(GY(1,3)*cos(t13) + BY(1,3)*sin(t13)) );
    h12 = V1*( -V1*BY(1,1) ...                                                 % Q1 (virtual)
             + V2*(GY(1,2)*sin(t12) - BY(1,2)*cos(t12)) ...
             + V3*(GY(1,3)*sin(t13) - BY(1,3)*cos(t13)) );


    zest = [h1; h2; h3; h4; h5; h6; h7; h8; h9; h10; h11; h12];   % zest = h(x)


    % ---- Jacobiana H(x) ----
    % --- Termos A e D de cada linha individual (G, B da linha) ---
    A21  = G12a*cos(t21) + B12a*sin(t21);   D21  = G12a*sin(t21) - B12a*cos(t21);  % linhas a e b (identicas)
    A31  = G13 *cos(t31) + B13 *sin(t31);   D31  = G13 *sin(t31) - B13 *cos(t31);
    A32  = G23 *cos(t32) + B23 *sin(t32);   D32  = G23 *sin(t32) - B23 *cos(t32);
    A12a = G12a*cos(t12) + B12a*sin(t12);   D12a = G12a*sin(t12) - B12a*cos(t12);


    % --- Termos A e D da Ybus completa (para as medidas virtuais) ---
    A12Y = GY(1,2)*cos(t12) + BY(1,2)*sin(t12);   D12Y = GY(1,2)*sin(t12) - BY(1,2)*cos(t12);
    A13Y = GY(1,3)*cos(t13) + BY(1,3)*sin(t13);   D13Y = GY(1,3)*sin(t13) - BY(1,3)*cos(t13);


    H = zeros(12, 5);


    % --- V1, V2, V3 ---
    H(1,3) = 1;
    H(2,4) = 1;
    H(3,5) = 1;


    % --- P21a  (i = 2, k = 1; t2 fixo) ---
    H(4,1) =  V2*V1*D21;
    H(4,2) =  0;
    H(4,3) =  V2*A21;
    H(4,4) =  V1*A21 - 2*G12a*V2;
    H(4,5) =  0;


    % --- P21b  (identica a P21a, com G12b, B12b) ---
    H(5,1) =  V2*V1*D21;
    H(5,2) =  0;
    H(5,3) =  V2*A21;
    H(5,4) =  V1*A21 - 2*G12b*V2;
    H(5,5) =  0;


    % --- P31  (i = 3, k = 1) ---
    H(6,1) =  V3*V1*D31;
    H(6,2) = -V3*V1*D31;
    H(6,3) =  V3*A31;
    H(6,4) =  0;
    H(6,5) =  V1*A31 - 2*G13*V3;


    % --- P32  (i = 3, k = 2; t2 fixo) ---
    H(7,1) =  0;
    H(7,2) = -V3*V2*D32;
    H(7,3) =  0;
    H(7,4) =  V3*A32;
    H(7,5) =  V2*A32 - 2*G23*V3;


    % --- Q12a  (i = 1, k = 2; t2 fixo) ---
    H(8,1) =  V1*V2*A12a;
    H(8,2) =  0;
    H(8,3) =  V2*D12a + 2*B12a*V1;
    H(8,4) =  V1*D12a;
    H(8,5) =  0;


    % --- Q32  (i = 3, k = 2; t2 fixo) ---
    H(9,1) =  0;
    H(9,2) =  V3*V2*A32;
    H(9,3) =  0;
    H(9,4) =  V3*D32;
    H(9,5) =  V2*D32 + 2*B23*V3;


    % --- Qgbat = (BY11 + BY12 + BY13)*V1^2 ---
    H(10,3) = 2*(BY(1,1) + BY(1,2) + BY(1,3))*V1;


    % --- P1 virtual: P1 = V1*[V1*GY11 + V2*A12Y + V3*A13Y] ---
    H(11,1) = -V1*V2*D12Y - V1*V3*D13Y;
    H(11,2) =  V1*V3*D13Y;
    H(11,3) =  2*V1*GY(1,1) + V2*A12Y + V3*A13Y;
    H(11,4) =  V1*A12Y;
    H(11,5) =  V1*A13Y;


    % --- Q1 virtual: Q1 = V1*[-V1*BY11 + V2*D12Y + V3*D13Y] ---
    H(12,1) =  V1*V2*A12Y + V1*V3*A13Y;
    H(12,2) = -V1*V3*A13Y;
    H(12,3) = -2*V1*BY(1,1) + V2*D12Y + V3*D13Y;
    H(12,4) =  V1*D12Y;
    H(12,5) =  V1*D13Y;




    % ---------- 2) restringir as medidas mantidas ----------
    zest_a = zest(ativas2);
    H_a    = H(ativas2, :);


    % ---------- 3) residuos, ganho, J e residuos normalizados ----------
    r2     = z_a - zest_a;
    G2     = H_a' * W2 * H_a;
    J2     = sum(r2.^2 ./ var_a);
    Omega2 = R2 - H_a * (G2 \ H_a');
    rN2    = r2 ./ sqrt(abs(diag(Omega2)));


    % ---------- 4) criterio de paragem: max|dx| < tol ----------
    if max(abs(dx2)) < tol
        break
    end
    if it2 >= maxit
        warning('Numero maximo de iteracoes (%d) atingido sem convergencia.', maxit);
        break
    end


    % ---------- 5) correcao: dx = [H'WH]^(-1) H'W [z - zest] ----------
    dx2 = G2 \ (H_a' * W2 * r2);
    x2  = x2 + dx2;
    it2 = it2 + 1;
    fprintf('(sem P21b) Iteracao %d:  max|dx| = %.4e\n', it2, max(abs(dx2)));
end
```

```matlabTextOutput
(sem P21b) Iteracao 1:  max|dx| = 2.1106e-01
(sem P21b) Iteracao 2:  max|dx| = 1.9119e-02
(sem P21b) Iteracao 3:  max|dx| = 3.5504e-04
(sem P21b) Iteracao 4:  max|dx| = 2.5970e-06
```

```matlab


% ==================================================================
%  RESULTADOS DA ULTIMA ITERACAO (sem P21b)
% ==================================================================
fprintf('\nConvergiu (sem P21b): max|dx| = %.4e < tol = %.1e  (%d iteracoes)\n\n', max(abs(dx2)), tol, it2);
```

```matlabTextOutput
Convergiu (sem P21b): max|dx| = 2.5970e-06 < tol = 1.0e-04  (4 iteracoes)
```

```matlab


idx_reais2 = find(~virtual2);
[rN_max2, q2] = max(abs(rN2(idx_reais2)));
nome_max2 = Names2{idx_reais2(q2)};


Resumo_vals2 = {rad2deg(x2(1)); rad2deg(x2(2)); x2(3); x2(4); x2(5); ...
                sprintf('%d / %d', m2, k2); it2; J2; chi2_crit2; ...
                sprintf('%.2f (%s)', rN_max2, nome_max2)};
disp_resumo2 = table(Resumo_vals2, 'RowNames', Resumo_nomes, 'VariableNames', {'Valor'})
```

| |Valor|
|:--:|:--:|
|1 theta1 (graus)|3.0812|
|2 theta3 (graus)|10.9980|
|3 V1 (pu)|1.0912|
|4 V2 (pu)|0.9785|
|5 V3 (pu)|1.0989|
|6 m / k = m - n|'11 / 6'|
|7 Iteracoes|4|
|8 J(x)|15.3038|
|9 chi2_{k;0,05}|12.5916|
|10 max \|rN\| (medidas reais)|'3.84 (Q32)'|

```matlab


if J2 >= chi2_crit2
    fprintf('Teste chi^2: J = %.4f >= %.4f  ->  HA dados errados (hipotese aceite a %d%%)\n\n', J2, chi2_crit2, round((1-alfa)*100));
else
    fprintf('Teste chi^2: J = %.4f < %.4f  ->  nao se pode afirmar que haja dados errados\n\n', J2, chi2_crit2);
end
```

```matlabTextOutput
Teste chi^2: J = 15.3038 >= 12.5916  ->  HA dados errados (hipotese aceite a 95%)
```

```matlab


% --- zest_i = h(x2) de TODAS as medidas (incluindo a removida) e rN_i das medidas utilizadas ---
rN2_full = NaN(numel(z), 1);   rN2_full(ativas2) = rN2;      % NaN = medida removida
disp_estimados2 = array2table([z, zest, z - zest, rN2_full, double(ativas2)], ...
    'RowNames', MeasNames, ...
    'VariableNames', {'z', 'z_estimado', 'residuo', 'rN', 'usada'})
```

| |z|z_estimado|residuo|rN|usada|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 V1|1.1000|1.0912|0.0088|0.7326|1|
|2 V2|1|0.9785|0.0215|1.7808|1|
|3 V3|1.0700|1.0989|-0.0289|-2.3626|1|
|4 P21a|-0.7900|-0.7913|0.0013|0.0417|1|
|5 P21b|0.7900|-0.7913|1.5813|NaN|0|
|6 P31|1.7100|1.6916|0.0184|0.7056|1|
|7 P32|2.3400|2.3556|-0.0156|-0.7963|1|
|8 Q12a|0.9700|1.1298|-0.1598|-2.4018|1|
|9 Q32|1.2200|1.1104|0.1096|3.8383|1|
|10 Qgbat|2.6600|2.6194|0.0406|0.6635|1|
|11 P1v|0|-7.5835e-05|7.5835e-05|0.5151|1|
|12 Q1v|0|-8.7137e-05|8.7137e-05|1.6051|1|

```matlab


% ==================================================================
%  3.a ESTIMACAO: DADOS INICIAIS, REMOVENDO P21b E Q32
% ==================================================================
ativas3 = true(numel(z), 1);
ativas3(strcmp(MeasNames, 'P21b')) = false;    % 1.a medida removida: P21b (erro de sinal, ver 2.a estimacao)
ativas3(strcmp(MeasNames, 'Q32'))  = false;    % 2.a medida removida: Q32
% (para remover SO a Q32, apagar a linha de P21b acima)
 
Names3   = MeasNames(ativas3);
z_a      = z(ativas3);
var_a    = var(ativas3);
virtual3 = virtual(ativas3);
W3 = diag(1 ./ var_a);                         % pesos das medidas mantidas
R3 = diag(var_a);
m3 = nnz(ativas3);                             % numero de medidas utilizadas
k3 = m3 - n;                                   % graus de liberdade
if exist('chi2inv', 'file')
    chi2_crit3 = chi2inv(1 - alfa, k3);
else
    chi2_crit3 = tabela_chi2(k3);
end
 
x3  = [0; 0; 1; 1; 1];                         % flat start
dx3 = Inf(n, 1);
it3 = 0;
 
while true
    % ---------- 1) zest = h(x) e H(x) no estado atual x3 (todas as 12 linhas) ----------
    t1 = x3(1);  t3 = x3(2);  V1 = x3(3);  V2 = x3(4);  V3 = x3(5);
 
    % --- Diferencas angulares ---
    t12 = t1 - t2;   t21 = t2 - t1;
    t31 = t3 - t1;   t13 = t1 - t3;
    t32 = t3 - t2;
 
    % --- Funcoes de medida ---
    h1  = V1;                                                                  % V1
    h2  = V2;                                                                  % V2
    h3  = V3;                                                                  % V3
    h4  = V2*V1*(G12a*cos(t21) + B12a*sin(t21)) - G12a*V2^2;                   % P21a
    h5  = V2*V1*(G12b*cos(t21) + B12b*sin(t21)) - G12b*V2^2;                   % P21b
    h6  = V3*V1*(G13 *cos(t31) + B13 *sin(t31)) - G13 *V3^2;                   % P31
    h7  = V3*V2*(G23 *cos(t32) + B23 *sin(t32)) - G23 *V3^2;                   % P32
    h8  = V1*V2*(G12a*sin(t12) - B12a*cos(t12)) + B12a*V1^2;                   % Q12a
    h9  = V3*V2*(G23 *sin(t32) - B23 *cos(t32)) + B23 *V3^2;                   % Q32
    h10 = (BY(1,1) + BY(1,2) + BY(1,3))*V1^2;                                  % Qgbat = Bcap*V1^2
    h11 = V1*( V1*GY(1,1) ...                                                  % P1 (virtual)
             + V2*(GY(1,2)*cos(t12) + BY(1,2)*sin(t12)) ...
             + V3*(GY(1,3)*cos(t13) + BY(1,3)*sin(t13)) );
    h12 = V1*( -V1*BY(1,1) ...                                                 % Q1 (virtual)
             + V2*(GY(1,2)*sin(t12) - BY(1,2)*cos(t12)) ...
             + V3*(GY(1,3)*sin(t13) - BY(1,3)*cos(t13)) );
 
    zest = [h1; h2; h3; h4; h5; h6; h7; h8; h9; h10; h11; h12];   % zest = h(x)
 
    % ---- Jacobiana H(x) ----
    % --- Termos A e D de cada linha individual (G, B da linha) ---
    A21  = G12a*cos(t21) + B12a*sin(t21);   D21  = G12a*sin(t21) - B12a*cos(t21);  % linhas a e b (identicas)
    A31  = G13 *cos(t31) + B13 *sin(t31);   D31  = G13 *sin(t31) - B13 *cos(t31);
    A32  = G23 *cos(t32) + B23 *sin(t32);   D32  = G23 *sin(t32) - B23 *cos(t32);
    A12a = G12a*cos(t12) + B12a*sin(t12);   D12a = G12a*sin(t12) - B12a*cos(t12);
 
    % --- Termos A e D da Ybus completa (para as medidas virtuais) ---
    A12Y = GY(1,2)*cos(t12) + BY(1,2)*sin(t12);   D12Y = GY(1,2)*sin(t12) - BY(1,2)*cos(t12);
    A13Y = GY(1,3)*cos(t13) + BY(1,3)*sin(t13);   D13Y = GY(1,3)*sin(t13) - BY(1,3)*cos(t13);
 
    H = zeros(12, 5);
 
    % --- V1, V2, V3 ---
    H(1,3) = 1;
    H(2,4) = 1;
    H(3,5) = 1;
 
    % --- P21a  (i = 2, k = 1; t2 fixo) ---
    H(4,1) =  V2*V1*D21;
    H(4,2) =  0;
    H(4,3) =  V2*A21;
    H(4,4) =  V1*A21 - 2*G12a*V2;
    H(4,5) =  0;
 
    % --- P21b  (identica a P21a, com G12b, B12b) ---
    H(5,1) =  V2*V1*D21;
    H(5,2) =  0;
    H(5,3) =  V2*A21;
    H(5,4) =  V1*A21 - 2*G12b*V2;
    H(5,5) =  0;
 
    % --- P31  (i = 3, k = 1) ---
    H(6,1) =  V3*V1*D31;
    H(6,2) = -V3*V1*D31;
    H(6,3) =  V3*A31;
    H(6,4) =  0;
    H(6,5) =  V1*A31 - 2*G13*V3;
 
    % --- P32  (i = 3, k = 2; t2 fixo) ---
    H(7,1) =  0;
    H(7,2) = -V3*V2*D32;
    H(7,3) =  0;
    H(7,4) =  V3*A32;
    H(7,5) =  V2*A32 - 2*G23*V3;
 
    % --- Q12a  (i = 1, k = 2; t2 fixo) ---
    H(8,1) =  V1*V2*A12a;
    H(8,2) =  0;
    H(8,3) =  V2*D12a + 2*B12a*V1;
    H(8,4) =  V1*D12a;
    H(8,5) =  0;
 
    % --- Q32  (i = 3, k = 2; t2 fixo) ---
    H(9,1) =  0;
    H(9,2) =  V3*V2*A32;
    H(9,3) =  0;
    H(9,4) =  V3*D32;
    H(9,5) =  V2*D32 + 2*B23*V3;
 
    % --- Qgbat = (BY11 + BY12 + BY13)*V1^2 ---
    H(10,3) = 2*(BY(1,1) + BY(1,2) + BY(1,3))*V1;
 
    % --- P1 virtual: P1 = V1*[V1*GY11 + V2*A12Y + V3*A13Y] ---
    H(11,1) = -V1*V2*D12Y - V1*V3*D13Y;
    H(11,2) =  V1*V3*D13Y;
    H(11,3) =  2*V1*GY(1,1) + V2*A12Y + V3*A13Y;
    H(11,4) =  V1*A12Y;
    H(11,5) =  V1*A13Y;
 
    % --- Q1 virtual: Q1 = V1*[-V1*BY11 + V2*D12Y + V3*D13Y] ---
    H(12,1) =  V1*V2*A12Y + V1*V3*A13Y;
    H(12,2) = -V1*V3*A13Y;
    H(12,3) = -2*V1*BY(1,1) + V2*D12Y + V3*D13Y;
    H(12,4) =  V1*D12Y;
    H(12,5) =  V1*D13Y;
 
 
 
    % ---------- 2) restringir as medidas mantidas ----------
    zest_a = zest(ativas3);
    H_a    = H(ativas3, :);
 
    % ---------- 3) residuos, ganho, J e residuos normalizados ----------
    r3     = z_a - zest_a;
    G3     = H_a' * W3 * H_a;
    J3     = sum(r3.^2 ./ var_a);
    Omega3 = R3 - H_a * (G3 \ H_a');
    rN3    = r3 ./ sqrt(abs(diag(Omega3)));
 
    % ---------- 4) criterio de paragem: max|dx| < tol ----------
    if max(abs(dx3)) < tol
        break
    end
    if it3 >= maxit
        warning('Numero maximo de iteracoes (%d) atingido sem convergencia.', maxit);
        break
    end
 
    % ---------- 5) correcao: dx = [H'WH]^(-1) H'W [z - zest] ----------
    dx3 = G3 \ (H_a' * W3 * r3);
    x3  = x3 + dx3;
    it3 = it3 + 1;
    fprintf('(sem P21b e Q32) Iteracao %d:  max|dx| = %.4e\n', it3, max(abs(dx3)));
end
```

```matlabTextOutput
(sem P21b e Q32) Iteracao 1:  max|dx| = 2.2641e-01
(sem P21b e Q32) Iteracao 2:  max|dx| = 2.4946e-02
(sem P21b e Q32) Iteracao 3:  max|dx| = 2.0813e-04
(sem P21b e Q32) Iteracao 4:  max|dx| = 1.0976e-06
```

```matlab
 
% ==================================================================
%  RESULTADOS DA ULTIMA ITERACAO (sem P21b e Q32)
% ==================================================================
fprintf('\nConvergiu (sem P21b e Q32): max|dx| = %.4e < tol = %.1e  (%d iteracoes)\n\n', max(abs(dx3)), tol, it3);
```

```matlabTextOutput
Convergiu (sem P21b e Q32): max|dx| = 1.0976e-06 < tol = 1.0e-04  (4 iteracoes)
```

```matlab
 
idx_reais3 = find(~virtual3);
[rN_max3, q3] = max(abs(rN3(idx_reais3)));
nome_max3 = Names3{idx_reais3(q3)};
 
Resumo_vals3 = {rad2deg(x3(1)); rad2deg(x3(2)); x3(3); x3(4); x3(5); ...
                sprintf('%d / %d', m3, k3); it3; J3; chi2_crit3; ...
                sprintf('%.2f (%s)', rN_max3, nome_max3)};
disp_resumo3 = table(Resumo_vals3, 'RowNames', Resumo_nomes, 'VariableNames', {'Valor'})
```

| |Valor|
|:--:|:--:|
|1 theta1 (graus)|3.1545|
|2 theta3 (graus)|11.5465|
|3 V1 (pu)|1.0980|
|4 V2 (pu)|1.0014|
|5 V3 (pu)|1.0735|
|6 m / k = m - n|'10 / 5'|
|7 Iteracoes|4|
|8 J(x)|0.5609|
|9 chi2_{k;0,05}|11.0705|
|10 max \|rN\| (medidas reais)|'0.65 (P31)'|

```matlab
 
if J3 >= chi2_crit3
    fprintf('Teste chi^2: J = %.4f >= %.4f  ->  HA dados errados (hipotese aceite a %d%%)\n\n', J3, chi2_crit3, round((1-alfa)*100));
else
    fprintf('Teste chi^2: J = %.4f < %.4f  ->  nao se pode afirmar que haja dados errados\n\n', J3, chi2_crit3);
end
```

```matlabTextOutput
Teste chi^2: J = 0.5609 < 11.0705  ->  nao se pode afirmar que haja dados errados
```

```matlab
 
% --- zest_i = h(x3) de TODAS as medidas (incluindo a removida) e rN_i das medidas utilizadas ---
rN3_full = NaN(numel(z), 1);   rN3_full(ativas3) = rN3;      % NaN = medida removida
disp_estimados3 = array2table([z, zest, z - zest, rN3_full, double(ativas3)], ...
    'RowNames', MeasNames, ...
    'VariableNames', {'z', 'z_estimado', 'residuo', 'rN', 'usada'})
```

| |z|z_estimado|residuo|rN|usada|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1 V1|1.1000|1.0980|0.0020|0.1660|1|
|2 V2|1|1.0014|-0.0014|-0.1330|1|
|3 V3|1.0700|1.0735|-0.0035|-0.3422|1|
|4 P21a|-0.7900|-0.7953|0.0053|0.1722|1|
|5 P21b|0.7900|-0.7953|1.5853|NaN|0|
|6 P31|1.7100|1.6929|0.0171|0.6531|1|
|7 P32|2.3400|2.3501|-0.0101|-0.5175|1|
|8 Q12a|0.9700|0.9565|0.0135|0.2755|1|
|9 Q32|1.2200|0.5613|0.6587|NaN|0|
|10 Qgbat|2.6600|2.6524|0.0076|0.1253|1|
|11 P1v|0|-7.4115e-05|7.4115e-05|0.5063|1|
|12 Q1v|0|2.4580e-05|-2.4580e-05|-0.5383|1|