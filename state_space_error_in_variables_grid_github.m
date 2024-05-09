%==========================================================================
% Companion code (2/2) for the paper 
% "Controller synthesis for input-state data with measurement errors."
% Evaluates the performance of the approaches with an energy and an
% instantaneous bound.
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
figIdx=0;                   % numbers figures in ascending order
nDataSets=20;
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

%% selection of parameters (\theta,T)
thetaVec  = [1e-6 sqrt(1e-11) 1e-5 sqrt(1e-9) 1e-4 sqrt(1e-7) 1e-3];
nThetaVec = length(thetaVec);
Tvec      = [20 40 60 80 100 120 140 160 180 200];
nTvec     = length(Tvec);
Tmax      = max(Tvec);

[THETAVEC_en, TVEC_en] = meshgrid(thetaVec,Tvec);
FEASPERC_en            = zeros(size(THETAVEC_en)); % initialization
[THETAVEC_in, TVEC_in] = meshgrid(thetaVec,Tvec);
FEASPERC_in            = zeros(size(THETAVEC_in)); % initialization
fprintf(1,'Progress:    0%%');

for j=1:nDataSets
    for jj=1:nThetaVec
    %% for each error bound generate the data corresponding to the variable Tmax
        thetaCurr=thetaVec(jj);
        % initial condition
        x0_exp=randn(n,1);
        % measured input
        um_exp=randn(m,Tmax);
        % error ex with |ex|^2 <= ex_bar = \theta/3
        ex_exp     = randn(n,Tmax+1);
        ex_exp     = ex_exp./sqrt(sum(ex_exp.^2,1));
        scalRadius = sqrt(thetaCurr/3)*nthroot(rand(1,Tmax+1),n);
        ex_exp     = ex_exp.*scalRadius;
        % error eu with |eu|^2 <= eu_bar = \theta/3
        eu_exp     = randn(m,Tmax);
        eu_exp     = eu_exp./sqrt(sum(eu_exp.^2,1));
        scalRadius = sqrt(thetaCurr/3)*nthroot(rand(1,Tmax),m);
        eu_exp     = eu_exp.*scalRadius;

        x_exp=zeros(n,Tmax+1);
        x_exp(:,1)=x0_exp;
        u_exp=um_exp-eu_exp;
        for k=idx(0):idx(Tmax-1)
            x_exp(:,k+1)=A_star*x_exp(:,k)+B_star*u_exp(:,k);
        end
        xm_exp=x_exp+ex_exp;

        numTresh=1e4;
        if max(vecnorm(x_exp))>numTresh
            error(['Norm of the state is above ' num2str(numTresh) '. Select shorter duration.'])
        end
        
        for jjj=1:nTvec
            %% for each length T, cut generated data
            currT=Tvec(jjj);
            currXm1 = x_exp(:,idx(1):idx(currT));
            currXm0 = x_exp(:,idx(0):idx(currT-1));
            currUm0 = u_exp(:,idx(0):idx(currT-1));
            currTheta   = currT*thetaCurr*eye(2*n+m);
            currTheta11 = currTheta(1:n,1:n);
            currTheta12 = currTheta(1:n,(n+1):(n+n+m));
            currTheta22 = currTheta((n+1):(n+n+m),(n+1):(n+n+m));            
            currAscript = [currXm0; currUm0]*[currXm0; currUm0]'-currTheta22;
            currBscript = -currXm1*[currXm0; currUm0]'+currTheta12;
            currCscript = (currXm1*currXm1.')-currTheta11;

            %% solve with energy-based, count how many times feasible over all datasets
            % if feasible, increment FEASENPERC
            [P_en, K_en, outcome_en]=design_data_based(currAscript,currBscript,currCscript);
            if outcome_en==0 % for YALMIP, 0 means that the problem can be certified to be feasible
                FEASPERC_en(jjj,jj)=FEASPERC_en(jjj,jj)+1;
                currAbsEigVec=abs(eig(A_star+B_star*K_en));
                if ~all(currAbsEigVec<.999)
                    warning('Some energy-based closed-loop eigenvalues of the true system have modulus >0.999.');
                end                
            end

            %% solve with instantaneous-based, count how many times feasible over all datasets
            % if feasible, increment FEASINPERC
            [Pdb_inst, Kdb_inst, outcome_inst]=design_data_based_inst(currUm0,currXm0,currXm1,thetaCurr);
            if outcome_inst==0
                FEASPERC_in(jjj,jj)=FEASPERC_in(jjj,jj)+1;
                currAbsEigVec=abs(eig(A_star+B_star*Kdb_inst));
                if ~all(currAbsEigVec<.999)
                    warning('Some instantaneous-based closed-loop eigenvalues of the true system have modulus >0.999.');
                end
            end

        end
        % displays progress
        progr=round(((j-1)*nThetaVec*sum(Tvec)+(jj-1)*sum(Tvec)+sum(Tvec(1:jjj)))/(nDataSets*nThetaVec*sum(Tvec))*100);
        fprintf(1,'\b\b\b\b%3d%%',progr);        
    end
end
fprintf('\n')

FEASPERC_en=FEASPERC_en/nDataSets;
FEASPERC_in=FEASPERC_in/nDataSets;

%% visualization

% adds fictitious 0's and a 1 for the sake of the colorbar
FEASPERC_en(end+1,end+1)=1;
FEASPERC_in(end+1,end+1)=1;
% makes the corresponding grid consistent in terms of dimensions
THETAVEC_en(:,end+1) = THETAVEC_en(:,end)+THETAVEC_en(:,end-1);
THETAVEC_en(end+1,:) = THETAVEC_en(end,:);
TVEC_en(end+1,:)     = TVEC_en(end,:)+TVEC_en(end-1,:);
TVEC_en(:,end+1)     = TVEC_en(:,end);
THETAVEC_in(:,end+1) = THETAVEC_in(:,end)+THETAVEC_in(:,end-1);
THETAVEC_in(end+1,:) = THETAVEC_in(end,:);
TVEC_in(end+1,:)     = TVEC_in(end,:)+TVEC_in(end-1,:);
TVEC_in(:,end+1)     = TVEC_in(:,end);

figure(1)
clf;
surf(log10(THETAVEC_en),TVEC_en,FEASPERC_en)
hold on; grid on; box on;
view([0 90])
colorbar
xlim(log10(thetaVec([1 end]))) % removes from view the fictitious 0's and 1
ylim(Tvec([1 end]))            % removes from view the fictitious 0's and 1
xlabel('log_{10}(\theta)')
ylabel('T')

figure(2)
clf; 
surf(log10(THETAVEC_in),TVEC_in,FEASPERC_in)
hold on; grid on; box on;
colorbar
view([0 90])
xlim(log10(thetaVec([1 end])))
ylim(Tvec([1 end]))

xlabel('log_{10}(\theta)')
ylabel('T')

%% functions
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
    opti    = sdpsettings('verbose',0);
    diagnos = optimize(constr,0,opti);

    outcome = diagnos.problem;
    Pv      = value(P);
    Yv      = value(Y);
    Kv      = Yv/Pv;    
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
    opti    = sdpsettings('verbose',0);
    diagnos = optimize(constr,0,opti);

    outcome = diagnos.problem;
    Pv      = value(P);
    Yv      = value(Y);
    Kv      = Yv/Pv;

end