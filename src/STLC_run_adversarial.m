function [Sys, params] = STLC_run_adversarial_mpt(Sys, controller)
% STLC_run_adversarial      runs a receding horizon control problem 
%                           for the system described by Sys, using the
%                           provided controller and adversary optimizer
%                           objects, using the horizon defined by Sys.L
%                           
% Input: 
%       Sys: an STLC_lti instance
%       controller: a YALMIP opmitizer object representing the system's 
%                   optimization problem
%       adversary: a YALMIP opmitizer object representing the adversarial  
%                  environment's optimization problem
%
% Output: 
%       Sys: modified with additional system_data
%       params: controller data
%
% :copyright: TBD
% :license: TBD

global StopRequest
StopRequest=0;

rob = Sys.min_rob;

%% Time
sys = Sys.sys;
ts=Sys.ts; % sampling time
L=Sys.L;  % horizon (# of steps)
time = Sys.time; % time for the date
time_d = (0:2*L-1)*ts; % discretized time for the controller

%% System dimensions and variables
Sys.sysd = c2d(sys,ts);

x0 = Sys.x0;
nu=Sys.nu;
nx=Sys.nx;
nw=Sys.nw;
ny=Sys.ny;

M = Sys.bigM; % big M


% get environment spec constraints (so we only compute them once)
[X_adv, Y_adv, U_adv, W_adv, Pstl, Fstl] = getSTL( Sys );

% reference disturbance
if isempty(Sys.Wref)
    Sys.Wref = 0*time;
end
Wref = Sys.Wref;
for iwx=1:nw
    Wrefn(iwx,:) = interp1( time , Wref(iwx,:)', time_d)';
end


%% Make data nb_stages times longer
nb_stages=Sys.nb_stages;
ntime = zeros(1, nb_stages*numel(time));
for istage = 0:nb_stages-1
    ntime(istage*numel(time)+1:(istage+1)*numel(time))= time+istage*(time(end)+time(2)) ;
end
time = ntime;
Wref = repmat(Wref,1,Sys.nb_stages);
Wn = Wrefn;


%% Initialize discrete data for the controller and environment

donen = zeros(1,2*L-1); % done(1:k) = 1 iff everything has been computed up to step k
pn = -1*M*ones(1,L);    % for activating robustness constraints
Un = zeros(nu,2*L-1);
Xn = zeros(max(nx,1),2*L);
if (nx>0)
    Xn(:,1) = x0;           % only X0 is already computed
end
pn(1) = rob;


Upred = zeros(nu,2*L-1);
Xpred = zeros(nx,2*L);

u_new = [];
w_new = [];
x_new = [];
y_new = [];
time_new = 0;
params = {};

% call solver for the first horizon
i_transient=1;
i_past = 1;

%% Init system and model data
Sys.system_data=struct;
Sys.model_data=struct;

Sys.system_data.time = [];
Sys.system_data.U = [];
Sys.system_data.X = x0;
Sys.system_data.Y = [];
Sys.system_data.W = [];

Sys.model_data.time = time_d;
Sys.model_data.X = repmat(0*time_d(1:end), [nx 1]);
Sys.model_data.Y = repmat(0*time_d(1:end-1), [ny 1]);
Sys.model_data.U = repmat(0*time_d(1:end-1), [nu 1]);


% compute initial environment using explicit MPC
[sol] =  compute_w();
if StopRequest==1
    return;
end
compute_input();


Sys.model_data.W = Wn;

u_new = Upred(:,i_transient);

w_new = rand_w();

[x_new, y_new] = system_step(Sys,x0, u_new, w_new);
params{end+1} = {time_d,donen,pn,Xn,Un,Wn,Wrefn};
i_past =  i_past+1;

update_hist_data();
Sys = update_plot(Sys);

%% loop
pause
if i_transient < L
    i_transient = i_transient+1;
end
while (time_d(end)+ts< time(end))
    % pause;
    x0 = x_new;
    time_new = time_new+ts;
    
    %% updates the model of controller and environment for the next horizon
    update_controller_data();
    
    %% compute input for the next horizon
    
    tic;
    [sol] = compute_w();
    toc;
    if StopRequest==1
        break;
    end
    tic;
    compute_input();
    toc;
    
    %% update states
    u_new = Upred(:,i_transient);
    w_new = rand_w();
    
    [x_new, y_new] =  system_step(Sys, x0, u_new, w_new);
    params{end+1} = {time_d,donen,pn,Xn,Un,Wn,Wrefn};
    i_past =  i_past+1;

    %% Update plots
    update_hist_data();
    Sys= update_plot(Sys);
    
    if i_transient < L
        i_transient = i_transient+1;
    end
    drawnow;
    %  pause;
    if StopRequest
        break;
    end
end

    function compute_input()
        
            Uopt = {};
            Xopt = {};
            Wopt = {};
            Jopt = [];
            
            % controller's turn
            for i=1:size(sol{1}.Pn,2)
                i
                size(sol{1}.Pn,2)
                A = full(sol{1}.Fi{i})
                B = full(sol{1}.Gi{i})
                [sol_control, errorflag1] = controller{{donen,pn,Xn,Un,A,B}};
                if(errorflag1==0)  % found a good control
                   disp(['Yalmip: ' yalmiperror(errorflag1)])
                   %disp(['Yalmip: ' 'Found a good control input'])
                   Uopt{i} = sol_control{1};
                   Xopt{i} = sol_control{2};
                   Wopt{i} = (sol{1}.Fi{i} * Uopt{i}' + sol{1}.Gi{i})';
                   Jopt = [Jopt;min(sol_control{3})];
                elseif (errorflag1==1 || errorflag1==15||errorflag1==12)  % some error, infeasibility or else
                   disp(['Yalmip error (disturbance too bad ?): ' yalmiperror(errorflag1)]); % probably there is no controller for this w
                   StopRequest=1;
                   return;
                else
                   disp(['Yalmip error: ' yalmiperror(errorflag1)]); % some other error
                   return;
                end
               
            end

            Upred = Uopt{find(Jopt==min(Jopt),1)};
            Xpred = Xopt{find(Jopt==min(Jopt),1)};
            %Wn = Wopt{find(Jopt==min(Jopt),1)};
                
    end

    function update_hist_data()
        
        Sys.system_data.time(end+1) = time_new;
        Sys.system_data.U(:,end+1) = u_new;
        Sys.system_data.X(:,end+1) = x_new;
        Sys.system_data.Y(:,end+1) = y_new;
        Sys.system_data.W(:,end+1) = w_new;
        
        Sys.model_data.time = time_d;
        Sys.model_data.X = double(Xpred);
        Sys.model_data.Y = [Sys.sysd.C*double(Xpred(:,1:end-1)) + Sys.sysd.D*[double(Upred); double(Wn(:,1:end-1))]];
        Sys.model_data.U = double(Upred);
        Sys.model_data.W = Wn;
        
    end

    function update_controller_data()
        
        if (i_past>= L+1)
            time_d = time_d+ts; % move forward one time step
        end
        
        % Note: the following assumes that system and model data are
        % both sampled with ts - should be changed eventually
        
        % reinitialize disturbance to the reference
        for wx=1:nw
            Wrefn(wx,:) = interp1( time , Wref(wx,:)', time_d)';
        end
        Wn = Wrefn;
        
        if i_transient<L
            donen(1:i_transient-1) = 1;  % we reuse what has been computed at the previous step
            Un(:,1:i_transient-1) = Sys.system_data.U(:,1:i_transient-1); %  previously computed inputs
            Wn(:,1:i_transient-1) = Sys.system_data.W(:,1:i_transient-1);
            Xn(:,1:i_transient) = Sys.system_data.X(:, 1:i_transient);
            pn(i_transient) = rob;
        else
            pn(1:L) = rob*ones(1,L);
            donen(1:L-1) = 1;
            Un(:,1:L-1) = Sys.system_data.U(:,end-L+2:end);   %    previously computed inputs
            Xn(:,1:L) = Sys.system_data.X(:,end-L+1:end);     %    previously computed temperatures
            Wn(:,1:L-1) = Sys.system_data.W(:,end-L+2:end);
        end
        
        
        
    end

    function w_new= rand_w() 
        % computes a random disturbance between w bounds 
        dw = Sys.w_ub- Sys.w_lb;    
        w_new = Wrefn(:,i_transient)+ dw'.*(2*rand(Sys.nw,1)-1)/3;        
    end


    function [sol] = compute_w() 
        fprintf('Starting MPT call...')
        [sol] = STLC_get_adversary_mpt(Sys, X_adv, Y_adv, U_adv, W_adv, Pstl, Fstl, donen, Xn, Un, Wn, Wrefn);
        fprintf('Done with mpt call...');
    end


    function [X, Y, U, W, Pstl, Fstl ] = getSTL( Sys )
        % variables
        X = sdpvar(nx, 2*L);
        Y = sdpvar(ny, 2*L-1);
        U = sdpvar(nu, 2*L-1);
        W = sdpvar(nw, 2*L);

        %% STL formula 
        Fstl=[];         
        Pstl = [];

        varStd = struct('X',X,'Y',Y,'U',U, 'W', W);

        if isstruct(Sys.var)
            %remove overlapping fields from std
            var = rmfield(varStd, intersect(fieldnames(Sys.var), fieldnames(varStd)));
            keys = [fieldnames(var); fieldnames(Sys.var)];
            var = cell2struct([struct2cell(varStd); struct2cell(Sys.var)], keys, 1);
        else
            var = varStd;
        end

        stl_list= STLC_parse_stl_labels(Sys);

        Pphi=sdpvar(1,1);
        for i = 1:numel(stl_list)
            phi = STLformula('phi', stl_list{i});
            [Fphi, Pphi] = STL2MILP_robust(phi, 2*L, ts, var, M);
            Fstl = [Fstl Fphi];
            Pstl = [Pstl,Pphi];
        end
    end






end

function Stop()
global StopRequest;
StopRequest = true;
end

