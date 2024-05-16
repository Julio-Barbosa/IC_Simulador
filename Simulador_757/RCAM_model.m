function [XDOT] = RCAM_model(X,U)

% ------------------- Vetores de Estado e Controle -----------------------
x1 = X(1); %u
x2 = X(2); %v
x3 = X(3); %w
x4 = X(4); %p
x5 = X(5); %q
x6 = X(6); %r
x7 = X(7); %phi
x8 = X(8); %theta
x9 = X(9); %psi 

u1 = U(1); %aileron
u2 = U(2); %profundor
u3 = U(3); %leme
u4 = U(4); %motor1
u5 = U(5); %motor2

% ------------------- Características da Aeronave -------------------------

m = 120000; % Massa total do veículo (kg)

cbar = 6.6; % Corda média aerodinâmica (m)
lt = 24.8; % Distânica do CA da cauda ao do corpo (m)
S = 260; % Área da asa (m²)
St = 64; % Área da cauda (m²)

Xcg = 0.23*cbar; % posição x do CG
Ycg = 0; % posição y do CG
Xcg = 0.10*cbar; % posição z do CG

Xac = 0.12*cbar; % posição x do CA
Yac = 0; % posição y do CA
Xac = 0; % posição z do CA

% ---- Motores ----

Xapt1 = 0; % posição x da força do motor 1
Yapt1 = -7.94; % posição y da força do motor 1
Zapt1 = -1.9; % posição z da força do motor 1

Xapt2 = 0; % posição x da força do motor 2
Yapt2 = 7.94; % posição y da força do motor 2
Zapt2 = -1.9; % posição z da força do motor 2

% ---- Outras constantes ----

rho = 1.225; % densidade do ar
g = 9.81; % aceleração gravitacional
depsda = 0.25; % mudança no ângulo de dawnwash
alpha_L0 = -11.5*pi/180; % ângulo de sustentação zero
n = 5.5; % inclinação da região linear da sustentação
a3 = -768.5; % coeficiente de alpha^3
a2 = 609.2; % coeficiente de alpha^2
a1 = -155.2; % coeficiente de alpha^1
a0 = 15.212; % coeficiente de alpha^0
alpha_switch = 14.5*(pi/180); % alpha onde sustentação não é mais linear

% --------------- Limites de Controle / Saturações -----------------------

u1min = -25*pi/180;
u1max = 25*pi/180;

u2min = -25*pi/180;
u2max = 10*pi/180;

u3min = -30*pi/180;
u3max = 30*pi/180;

u4min = 0.5*pi/180;
u4max = 10*pi/180;

u5min = 0.5*pi/180;
u5max = 10*pi/180;

if u1>u1max
    u1 = u1max;
elseif u1<u1min
    u1 = u1min;
end

if u2>u2max
    u2 = u2max;
elseif u2<u2min
    u2 = u2min;
end

if u3>u3max
    u3 = u3max;
elseif u3<u3min
    u3 = u3min;
end

if u4>u4max
    u4 = u4max;
elseif u4<u4min
    u4 = u4min;
end

if u5>u5max
    u5 = u5max;
elseif u5<u5min
    u5 = u5min;
end

% --------------------- Variáveis Intermediárias -------------------------

% Velocidade do ar
Va = sqrt(x1^2 + x2^2 + x3^2);

% Calculo alpha e beta
alpha = atan2(x2,x1);
beta = asin(x2/Va);

% Pressão Dinâmica
Q = 0.5*rho*Va^2;

% Alguns vetores
wbe_b = [x4;x5;x6];
V_b = [x1;x2;x3];

% ---------------------------- Aerdinâmica -------------------------------

% Cálculo CL_wb
if alpha <= alpha_switch
    CL_wb = n*(alpha - alpha_L0);
else 
    CL_wb = a3*alpha^3 + a2*alpha^2 + a1*alpha^1 + a0;
end 

% Cálculo CL_t
epsilon = depsda*(alpha - alpha_L0);
alpha_t = alpha - epsilon + u2 + 1.3*x5*lt/Va;
CL_t = 3.1*(St/S)*alpha_t;

% Sustetação Total
CL = CL_wb + CL_t;

% Arrasto Total
CD = 0.13 + 0.07*(5.5*alpha + 0.654)^2;

% Força Lateral
CY = -1.6*beta + 0.24*u3;

% ---- Dimensionamento ----

% Cáculo das Forças
FA_s = [-CD*Q*S;
         CY*Q*S;
        -CL*Q*S];

% Rotacionando forças eixo corpo
C_bs = [cos(alpha) 0 -sin(alpha);
        0 1 0;
        sin(alpha) 0 cos(alpha)];
FA_b = C_bs*FA_s;

% ---- Momento Aerodinâmico ----

eta11 = -1.4*beta;
eta21 = -0.59 - (3.1*(St*lt)/(S*cbar))*(alpha - epsilon);
eta31 = (1 - alpha*(180/(15*pi)))*beta;

eta = [eta11;
       eta21;
       eta31];

dCMdx = (cbar/Va)*[-11 0 5;
                    0 (-4.03*(St*lt^2)/(S*cbar^2)) 0;
                    1.7 0 -11.5];

dCMdu = [-0.6 0 0.22;
          0 (-3.1*(St*lt)/(S*cbar)) 0;
          0 0 -0.63];

CMac_b = eta + dCMdx*wbe_b + dCMdu*[u1;u2;u3];

% Momento em torno do CA
MAac_b = CMac_b*Q*S*cbar;

% Momento em torno do CG
rcg_b = [Xcg;Ycg;Zcg];
rac_b = [Xac;Yac;Zac];
MAac_b = MAac_b + cross(FA_b,rcg_b - rac_b);

% ------------------------- Desempenho ----------------------------------- 

% Tração
F1 = u4*m*g;
F2 = u5*m*g;

FE1_B = [F1;0;0];
FE2_B = [F2;0;0];
FE_B = FE1_B + FE2_B;

% Momento
mew1 = [Xcg - Xapt1;
        Yapt1 - Ycg;
        Zcg - Zapt1];

mew2 = [Xcg - Xapt2;
        Yapt2 - Ycg;
        Zcg - Zapt2];

MEcg1_b = cross(mew1,FE1_b);
MEcg2_b = cross(mew2,FE2_b);
MEcg_b = MEcg1_b+MEcg2_b;

% ---- Efeitos Gravitacionais ----
g_b = [-g*sin(x8);
        g*cos(x8)*sin(x7);
        g*cos(x8)*cos(x7)];

Fg_b = m*g_b;

% ----------------------- Derivadas de Estado-----------------------------

% Matriz Inercial
Ib = m*[40.07 0 -2.0923;
        0 64 0;
        -2.0923 0 99.92];

% Iverso da Matriz Inercial
invIb = (1/m)*[0.0249836 0 0.000523151;
               0 0.015625 0;
               0.000523151 0 0.010019]; % não chamou a função inversa pra economizar tempo computacional

% Forças e Momentos em F_b
F_b = Fg_b + FE_b + FA_b;
x1to3dot = (1/m)*F_b - cross(wbe_b,V_b);

Mcg_b = MAcg_b + MEcg_b;
x4to6dot = invIb*(Mcg_b - cross(wbe_b,Ib*wbe_b));

% Cáculo ângulos Euler
H_phi = [1 sin(x7)*tan(x8) cos(x7)*tan(x8);
         0 cos(x7) -sin(x7);
         0 sin(x7)/cos(x8) cos(x7)/cos(x8)];

x7to9dot = H_phi*wbe_b;

% XDOT
XDOT = [x1to3dot
        x4to6dot
        x7to9dot]; 


