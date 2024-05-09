%==========================================================================
% Companion code (1/2) for the paper 
% "Controller synthesis for input-state data with measurement errors."
% For a dataset, solves the feasibility programs in Theorems 1 and 2.
% Author: Andrea Bisoffi
%==========================================================================
clearvars
clc

%% information for the user
fprintf(['This is the companion code for\n' ...
    '\t Bisoffi, A., Li, L., De Persis, C., & Monshizadeh, N. (2024).\n' ...
    '\t "Controller synthesis for input-state data with measurement errors."\n' ...
    '\t arXiv preprint arXiv:2402.04157.\n'])  
fprintf(['This code needs the toolbox YALMIP and has been run with solver MOSEK.\n' ...
    'YALMIP can be freely downloaded from https://yalmip.github.io/.\n' ...
    'MOSEK can be obtained from https://www.mosek.com/.\n'...
    'Other solvers can be used instead of MOSEK, but were not tested.\n'])

%% high-level parameters
eu_bar=5e-5;               % instantaneous bound on error |eu|^2<=eu_bar
ex_bar=5e-5;               % instantaneous bound on error |ex|^2<=ex_bar
T=200;                     % number of data
figIdx=0;                  % numbers the figures in ascending order
rng('default')
idx=@(x)(x+1);

%% system
% matrices A, B corresponding to a simple distillation column from p. 95 of
% P. Albertos, A. Sala, "Multivariable control systems: an engineering
% approach". Springer, 2004.
A_star=[0       0       0       0       0       0       0;
        1       0       0       0       0       0       0;
        0       0.0548  0.9217  0       0       0       0;
        0.4882  0       0       0.9024  0       0       0;
        0       0       0       0       0.8607  0       0;
        0       0       0       0       0       0.8607  0;
        0       0       0       0       0       0       0.9024];
B_star=[1       0       0;
        0       0       0;
        0       0       0;
        0       0       0;
        0.0348  0       0.2298;
        0       0.0808  0;
        0       0.3564  0];
[n, m]=size(B_star);

%% data collection experiment
% initial condition
x0_exp=randn(n,1);
% measured input
um_exp=randn(m,T);
% error ex with |ex|^2 <= ex_bar
ex_exp     = randn(n,T+1);
ex_exp     = ex_exp./sqrt(sum(ex_exp.^2,1));
scalRadius = sqrt(ex_bar)*nthroot(rand(1,T+1),n);
ex_exp     = ex_exp.*scalRadius;
% error eu with |eu|^2 <= eu_bar
eu_exp     = randn(m,T);
eu_exp     = eu_exp./sqrt(sum(eu_exp.^2,1));
scalRadius = sqrt(eu_bar)*nthroot(rand(1,T),m);
eu_exp     = eu_exp.*scalRadius;

t_exp=0:(T-1);
x_exp=zeros(n,T+1);

x_exp(:,1)=x0_exp;
u_exp=um_exp-eu_exp;
for k=idx(0):idx(T-1)
    x_exp(:,k+1)=A_star*x_exp(:,k)+B_star*u_exp(:,k);
end
xm_exp=x_exp+ex_exp;

numTresh=1e4;
if max(vecnorm(x_exp))>numTresh
    error(['Norm of the state is above ' num2str(numTresh) '. Select shorter duration.'])
end

X0      = x_exp(:,idx(0):idx(T-1));
U0      = u_exp(:,idx(0):idx(T-1));
Xm0     = xm_exp(:,idx(0):idx(T-1));
Xm1     = xm_exp(:,idx(1):idx(T));
Um0     = um_exp(:,idx(0):idx(T-1));
Theta   = T*(2*ex_bar+eu_bar)*eye(2*n+m);
Theta11 = Theta(1:n,1:n);
Theta12 = Theta(1:n,(n+1):(n+n+m));
Theta22 = Theta((n+1):(n+n+m),(n+1):(n+n+m));

figIdx=figIdx+1; figure(figIdx); clf; hold on; grid on; box on;
ph=plot(t_exp,Um0,'.-');
set(gca,'ColorOrderIndex',1)
plot(t_exp,U0,'.:')
title('\textbf{Data collection}: input','Interpreter','latex')
xlabel('Time','Interpreter','latex')
ylabel('$u^\mathrm{m}$, $u$','Interpreter','latex')
ah=gca; ah.TickLabelInterpreter='latex';
xlim(t_exp([1 end]))

figIdx=figIdx+1; figure(figIdx); clf; hold on; grid on; box on;
plot(t_exp,Xm0,'.-')
set(gca,'ColorOrderIndex',1)
plot(t_exp,X0,'.:')
title('\textbf{Data collection}: state','Interpreter','latex')
xlabel('Time','Interpreter','latex')
ylabel('$x^\mathrm{m}$, $x$','Interpreter','latex')
ah=gca; ah.TickLabelInterpreter='latex';
xlim(t_exp([1 end]))

%% solving the model-based case
[Pmb, Kmb]=design_model_based(A_star,B_star);

%% solving the data-based case
Ascript=[Xm0; Um0]*[Xm0; Um0]'-Theta22;
Bscript=-Xm1*[Xm0; Um0]'+Theta12;
Cscript=(Xm1*Xm1.')-Theta11;
Zscript=-(Bscript/Ascript);
Qscript=(Bscript/Ascript)*Bscript'-Cscript;
[Pdb, Kdb, outcome]=design_data_based(Ascript,Bscript,Cscript);
Kdb

%% solving the data-based case: instantaneous bound
theta=2*ex_bar+eu_bar;
[Pdb_inst, Kdb_inst, outcome_inst]=design_data_based_inst(Um0,Xm0,Xm1,theta);
Kdb_inst

%% illustrate closed-loop simulations
figIdx=figIdx+1; figure(figIdx); clf; hold on; grid on; box on;

% selects which design to simulate: mb is "model-based", db is "data-based"
% and dbi is "data-based instantaneous"
% flag='mb';
% flag='db';
flag='dbi';
switch flag
    case 'mb'
        K=Kmb;      
    case 'db'
        K=Kdb;
    case 'dbi'
        K=Kdb_inst;
end

T_cl=40;
t_cl=0:T_cl;
x_cl=zeros(n,T_cl);
x_cl(:,1)=x0_exp;

for k=idx(0):idx(T_cl-1)
    x_cl(:,k+1)=(A_star+B_star*K)*x_cl(:,k);
end

plot(t_cl,x_cl,'.-')
title('\textbf{Closed-loop validation}','Interpreter','latex')
xlabel('Time','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
ah=gca; ah.TickLabelInterpreter='latex';


%% functions
function [Pv, Kv, outcome, Yv]=design_model_based(A_star,B_star)
    numStrictZero=1e-3;

    [n, m]=size(B_star);
    Y     =sdpvar(m,n,'full');
    P     =sdpvar(n,n);
    
    lyapDecr      = blkvar;
    lyapDecr(1,1) = -P+numStrictZero*eye(n);
    lyapDecr(1,2) = -[A_star, B_star]*[P; Y];
    lyapDecr(2,2) = -P;
    lyapDecr      = sdpvar(lyapDecr);

    constr  = [P>=numStrictZero*eye(n); lyapDecr<=0];
    optio   = sdpsettings('verbose',0);
    diagnos = optimize(constr,0,optio);
    
    outcome = diagnos.problem;
    Pv      = value(P);
    Yv      = value(Y);
    Kv      = Yv/Pv;
    fprintf('Model-based design: %s.\n',diagnos.info)
end

function [Pv, Kv, outcome, Yv]=design_data_based(Ascript,Bscript,Cscript)
    numStrictZero=1e-4;

    n=size(Bscript,1);
    m=size(Ascript,1)-n;
    Y=sdpvar(m,n,'full');
    P=sdpvar(n,n);
    
    lyapDecr      = blkvar;
    lyapDecr(1,1) = -P-Cscript+numStrictZero*eye(n);
    lyapDecr(1,2) = 0;
    lyapDecr(1,3) = Bscript;
    lyapDecr(2,2) = -P;
    lyapDecr(2,3) = [P; Y]';
    lyapDecr(3,3) = -Ascript;
    lyapDecr      = sdpvar(lyapDecr);

    constr  = [P>=numStrictZero*eye(n); lyapDecr<=0];
    optio   = sdpsettings('verbose',0);
    diagnos = optimize(constr,0,optio);

    outcome = diagnos.problem;
    Pv      = value(P);
    Yv      = value(Y);
    Kv      = Yv/Pv;
    fprintf('Data-based design, energy bound: %s.\n',diagnos.info)
end

function [Pv, Kv, outcome, Yv]=design_data_based_inst(Um0,Xm0,Xm1,theta)
    numStrictZero=1e-4;

    [n, T] = size(Xm0);
    m      = size(Um0,1);
    P      = sdpvar(n,n);
    Y      = sdpvar(m,n,'full');
    tauVec = sdpvar(1,T,'full');
    
    dataMtx=zeros(n+n+m+n,n+n+m+n);
    for k=1:T
        singleDatum=[ Xm1(:,k);
                      -Xm0(:,k);
                      -Um0(:,k);
                      zeros(n,1)];
        sumTerm=(singleDatum*singleDatum.')-theta*blkdiag(eye(n+n+m),zeros(n,n));
        dataMtx=dataMtx+tauVec(k)*sumTerm;
    end
    lyapDecr=...
        [-P+numStrictZero*eye(n),   zeros(n,n), zeros(n,m),	zeros(n,n);
        zeros(n,n),                 P,          Y.',        zeros(n,n);
        zeros(m,n),                 Y,          zeros(m,m), Y;
        zeros(n,n),                 zeros(n,n), Y.',        -P]...
        - dataMtx;
    
    constr  = [P>=numStrictZero*eye(n), lyapDecr<=0, tauVec>=0];
    optio   = sdpsettings('verbose',0);
    diagnos = optimize(constr,0,optio);

    outcome = diagnos.problem;
    Pv      = value(P);
    Yv      = value(Y);
    Kv      = Yv/Pv;
    fprintf('Data-based design, instantaneous bound: %s.\n',diagnos.info)
end