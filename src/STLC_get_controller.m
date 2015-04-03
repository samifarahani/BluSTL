function controller = STLC_get_controller(Sys,enc)
% STLC_get_controller constructs the controller object for an STLC_lti instance
%
% Input:
% Sys: an STLC_lti instance
%
% Output:
% controller: a YALMIP optimizer object that solves the STL-constrained
% optimal control problem for Sys
%
% :copyright: TBD
% :license: TBD
if nargin < 2
enc='robust';
end
%% Time
ts=Sys.ts; % sampling time
L=Sys.L; % horizon (# of steps)
%% System dimensions and variables
nu=Sys.nu;
nx=Sys.nx;
nw=Sys.nw;
ny=Sys.ny;
% variables
X = sdpvar(nx, 2*L);
U = sdpvar(nu, 2*L-1);
Y = sdpvar(ny, 2*L-1);
% parameters
W = sdpvar(nw, 2*L);
done = binvar(1,2*L-1);
p = sdpvar(1,L);
Udone = sdpvar(nu,2*L-1);
Xdone = sdpvar(nx, 2*L);
Wref = sdpvar(nw, 2*L);

x_prev = sdpvar(nx,1);
% %% STL formula
% % Fstl = [];
% % varStd = struct('X',X,'Y', Y,'U',U, 'W', W);
% % if isstruct(Sys.var)
% % %remove overlapping fields from std
% % var = rmfield(varStd, intersect(fieldnames(Sys.var), fieldnames(varStd)));
% % keys = [fieldnames(var); fieldnames(Sys.var)];
% % var = cell2struct([struct2cell(varStd); struct2cell(Sys.var)], keys, 1);
% % else
% % var = varStd;
% % end
% % stl_list= STLC_parse_stl_labels(Sys);
M = Sys.bigM;
% % Pphi=sdpvar(1,1);
% % for i = 1:numel(stl_list)
% % phi = STLformula('phi', stl_list{i});
% % switch enc
% % case 'boolean'
% % [Fphi, Pphi] = STL2MILP_boolean(phi, 2*L, ts, var,M);
% % case 'robust'
% % [Fphi, Pphi] = STL2MILP_robust(phi, 2*L, ts, var,M);
% % case 'interval'
% % [Fphi, Pphi1, Pphi2] = STL2MILP_robust_interval(phi, 2*L, ts, var,M);
% % end
% % Fstl = [Fstl Fphi];
% % for j = 1:min(L, size(Pphi,2))
% % switch enc
% % case 'boolean'
% % Fstl = [Fstl Pphi(:,j) == 1]; % TODO this is specific to alw (phi), whatabout ev, until...
% % case 'robust'
% % Fstl = [Fstl Pphi(:,j)>= p(j)]; % TODO this is specific to alw (phi), whatabout ev, until...
% % case 'interval'
% % Fstl = [Fstl Pphi2(:,j)>= p(j)]; % TODO this is specific to alw (phi), whatabout ev, until...
% % for k=1:2*L
% % if k==1
% % Fstl = [Fstl, Pphi1(:,1)==Pphi2(:,2)];
% % else
% % % done values (history)
% % % if k is past (done(k)==1), upper and lower bounds are equal
% % Fstl = [Fstl, Pphi1(:,k) - (1-done(k-1))*M <= Pphi2(:,k) <= Pphi1(:,k) + (1-done(k-1))*M];
% % end
% % end
% % end
% % end
% % end
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
%% Dynamics constraints
Fdyn = [];
[Ad,Bd,Cd,Dd]=ssdata(Sys.sysd);
Bdu=Bd(:,1:nu);
Bdw=Bd(:,nu+1:end);
Ddu=Dd(:,1:nu);
Ddw=Dd(:,nu+1:end);
% Constraints for states (if any)
for k=1:2*L
if k==1
Fdyn = [Fdyn, X(:,1)==Xdone(:,1)];
else
% done values (history)
% if k is past (done(k)==1), use values in Tdone, otherwise use linear update
Fdyn = [Fdyn, Xdone(:,k) - (1-done(k-1))*M <= X(:,k) <= Xdone(:, k)+ (1-done(k-1))*M];
% not done values
Fdyn = [Fdyn, ((Ad*X(:,k-1) + Bdu*U(:,k-1) + Bdw*W( :, k-1 )) - done(k-1)*M) <= X(:,k) <= ((Ad*X(:,k-1) + Bdu*U(:,k-1) + Bdw*W( :, k-1 )) + done(k-1)*M)];
end
end
% Constraints for inputs
for k=1:2*L
if k>1
Fdyn = [Fdyn, Udone(:,k-1) - (1-done(k-1))*M <= U(:,k-1) <= Udone(:,k-1) + (1-done(k-1))*M];
end
end
% Constraints for outputs (if any)
for k=1:2*L-1
Fdyn = [Fdyn, Y(:,k) == Cd*X(:,k)+ Ddu*U(:,k) + Ddw *W(:,k)];
end
%% STL formula
wlb = Sys.w_lb;
wub = Sys.w_ub;

[Matrices,z,lbz,ubz] = get_MILP_matrices(Ad,Bdu,Bdw,nx,nu,nw,Wref,wlb,wub,L,M);
%% Objective function
% % switch enc
% % case 'boolean'
% % obj = get_objective(Sys,X,Y,U,W);
% % case 'robust'
% % obj = get_objective(Sys,X,Y,U,W, Pphi(:,L), Sys.lambda_rho);
% % case 'interval'
% % obj = get_objective(Sys,X,Y,U,W, Pphi1(:,L), Sys.lambda_rho);
% % end
% Solving the optimization as a MILP using lagrangian multiplier
    % Cmpc_i x(k-1) + Dmpc_i \tilde V(k) + Rmpc_i \tilde w(k) - Empc_i = beta_i + Lambda_i (S \tilde e - q) i=1,...,m
    % beta <= 0 and Lambda >= 0 (Lambda is a mtrix of size m*nq)
    % min_{u,Lambda,beta} J(\tilde u(k))= CC * \tilde u(k)
    % s.t. S^T * Lambda_i = Rmpc_i
    %      Cmpc_i x(k-1) + Dmpc_i \tilde V(k) - Empc_i = beta_i - Lambda_i q
    %      Hmpc \tilde V = Fmpc
    %      Lambda >= 0 , beta <= 0
    
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
    obj1 = CC2'*z;
    Fnew = [];
    klp = 0;
    for ilp=1:size(Dmpc1,1)
       Fnew = Fnew + set(Dmpc1(ilp,:)*U' + Dmpc2(ilp,:) * z - Gmpc(ilp) == beta(ilp) ...
                        - Lambda(klp+1:klp+length(qmpc))'*qmpc) + set(Smpc'*Lambda((klp+1:klp+length(qmpc))) ...
                        == Rmpc(ilp,:)') + set(Hmpc * V == Fmpc); 
       klp = klp + length(qmpc);
    end
    Fnew =Fnew + set(Lambda >=0) + set(beta <=0) + set(lbz <= z <= ubz);
%     sollp = solvesdp([Fnew, Fdyn, Fstl, Fuf],obj1);
%     x_fin = double(U);

options = Sys.solver_options;
param_controller = {done, Xdone, Udone, x_prev, Wref};
output_controller = {U,X};
controller = optimizer([Fnew, Fdyn, Fu],obj1,options,param_controller, output_controller);