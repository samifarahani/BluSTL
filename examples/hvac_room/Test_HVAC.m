% HVAC dual optimization
clear all;
clc;

Sys = hvac_room;
Sys = init_control(Sys,5,0, 600);

starttime = clock;
%% Time
sys = Sys.sys;
ts=Sys.ts; % sampling time
L=Sys.L;  % horizon (# of steps)
time = Sys.time; % time for the date
time_d = (0:2*L-1)*ts; % discretized time for the controller
%% System dimensions and variables
Sys.sysd = c2d(sys,ts);
nu=Sys.nu;
nx=Sys.nx;
nw=Sys.nw;
ny=Sys.ny;
M = Sys.bigM;
nz = 3;
kend = (time(end)+1)/ts;
x = zeros(nx,kend+1);        % state variable
u = zeros(nu,kend);      % control varibale
z_stl = zeros(nz,kend);
%e = 2*rand(ne,kend) - 1;      % uncertain vector e
x_init = zeros(nx,1);%x(:,1);           % initial state
x(:,1) = x_init;
x_prev = x_init;           % state at the time step k-1

ts=Sys.ts; % sampling time
L=Sys.L; % horizon (# of steps)
Wref = Sys.Wref;
Np = L;
Sys.model_data.time = time_d;

U = sdpvar(1,nu*(2*L-1));

%% dynamincs
[Ad,Bd,Cd,Dd]=ssdata(Sys.sysd);
Bdu=Bd(:,1:nu);
Bdw=Bd(:,nu+1:end);
Ddu=Dd(:,1:nu);
Ddw=Dd(:,nu+1:end);

%% Input constraints
Fu = [];
% Bounds
for iu = 1:nu
Fu = [ Fu, Sys.u_lb(iu) <= U(iu,:) <= Sys.u_ub(iu)] ; % bounds constraints on u
end
% Bounded variability
delta_not_inf = 0;
for iu = 1:nu
dif = sdpvar(nu,2*L-2);
F_dif = dif(:,1:2*L-2) == U(:,2:2*L-1) - U(:,1:2*L-2);
if (Sys.u_delta(iu) < Sys.u_ub(iu)- Sys.u_lb(iu))
Fdif = [F_dif, -Sys.u_delta(iu) <= dif <= Sys.u_delta(iu)];
delta_not_inf = 1;
end
end
if delta_not_inf
Fu = [Fu Fdif];
end

%% STL formula
wlb = Sys.w_lb;
wub = Sys.w_ub;
k = 1;    % time step counter 
 for k=1:kend  
     k
     for wx=1:nw
        Wrefn(wx,:) = interp1( time , Wref(wx,:)', time_d)';
     end
    [Matrices,z,lbz,ubz] = get_MILP_matrices(Ad,Bdu,Bdw,nx,nu,nw,Wrefn,wlb,wub,L,M);

    Cmpc = Matrices.Cmpc;
    Dmpc1 = Matrices.Dmpc1(:,1:end-1);
    Dmpc2 = Matrices.Dmpc2;
    Rmpc = Matrices.Rmpc;
    Empc = Matrices.Empc;
    Hmpc = Matrices.Hmpc;
    Fmpc = Matrices.Fmpc;
    Smpc = Matrices.Smpc;
    qmpc = Matrices.qmpc;
    CC2 = Matrices.CC2;
    Gmpc = Empc - Cmpc * x_prev;
    V = [U';z];
    Lambda1 = sdpvar(size(Dmpc1,1),length(qmpc));
    Lambda = reshape(Lambda1.',[],1);
    beta = sdpvar(size(Dmpc1,1),1);
%     Lambda = sdpvar(size(Dmpc1,1),1);
    obj1 = CC2'*z;
%     obj1 = CC1 * U + 
    Fnew = [];
    klp = 0;
    for ilp=1:size(Dmpc1,1)
       Fnew = Fnew + set(Dmpc1(ilp,:)*U' + Dmpc2(ilp,:) * z - Gmpc(ilp) == beta(ilp) ...
                        - Lambda(klp+1:klp+length(qmpc))'*qmpc) + set(Smpc'*Lambda((klp+1:klp+length(qmpc))) ...
                        == Rmpc(ilp,:)') + set(Hmpc * z == Fmpc); 
       klp = klp + length(qmpc);
    end
    Fnew =Fnew + set(Lambda >=0) + set(beta <=0) + set(lbz <= z <= ubz);
    sollp = solvesdp([Fnew,Fu],obj1);
    x_fin = double(U);
    z_fin = double(z);
    
    u(:,k) = x_fin(1:nu);
    z_stl(:,k) = z_fin(1:nz);
    % computes a random disturbance between w bounds 
    dw = Sys.w_ub- Sys.w_lb;    
    w_new = Wrefn(:,1)+ dw'.*(2*rand(Sys.nw,1)-1)/3;        
    %JJ(k) = CC(1:nu) * u(:,k) + CC(Nc*nu+1:Nc*nu+nd)* delta(:,k) + CC(Nc*nu+Np*nd+1:Nc*nu+Np*nd+nz)*z(:,k);
    x(:,k+1) = Ad*x_prev + Bdu*u(:,k) + Bdw*w_new;
    % x_next = A*x_prev + B1*u(:,k+1) + B2*delta(:,k+1) + B3*z(:,k+1) + B4*e(:,k+1);  % for uncertain system
    
    x_prev = x(:,k+1);
 end
 endtime=clock;
totaltime = etime(endtime,starttime)