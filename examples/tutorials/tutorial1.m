%% Defining the plant dynamics 
% The toolbox is organized around one main class, called STLClti. An
% STLClti object is primarily a continuous Linear Time Invariant (LTI) 
% system with inputs, outputs and disturbances. Hence, a constructors for this 
% class takes matrices to define such an LTI. We first define and A and B
% matrices for state evolution:
A = [0 1 ;
     0 0];
Bu = [0;1];

%%
% Later on, we will use a disturbance signals so we need to define a Bw 
% matrice. This signal will not influence the state dynamics, though, so we
% set Bw to be 0. 
Bw = [0;0]; 

%%
% Next we define the output dynamics, i.e., C, Du and Dw matrices. Here we
% have a single output $y(t) = x1(t)+w(t)$.
C = [1 0];
Du = 0;
Dw = 0;

%%
% Now we can call the main constructor of STLC_lti class. 
Sys= STLC_lti(A,Bu,Bw,C,Du,Dw); 

%%
% In the next section, we will define the different settings for the 
% control synthesis experiment. Before that, we define some initial state:
Sys.x0= [1 ; 1];

%% Defining the controller
% We start by defining the time instants for the whole experiment, the discrete time
% step ts for the controller and the horizon L in number of time steps.  
Sys.time = 0:.1:10; 
Sys.ts=.2; % sampling time for controller
Sys.L=10;  % horizon is 2s in that case

%%
% Next we declare some constraints on control inputs, here, lower and upper
% bounds:
Sys.u_ub = 10;  % upper bound on u 
Sys.u_lb = -10; % lower bound on u

%%
% Then the following define a signal temporal logic formula to be satisfied
% by the system. Note that times in the temporal operators are continuous,
% not discrete steps. 
Sys.stl_list = {'ev_[0,1.] alw_[0,0.5] ( abs(y1(t)-w1(t)) < 0.1)'};

%%
% Now we are ready to compile the controller for our problem. 
controller = get_controller(Sys)

%%
% Note that by default, the objective function will minimize the 1-norm of
% the input. 

%% Testing the controller
%  The simplest mode to run our system with the newly created controller is in
%  open loop. This is done with the following command:
run_open_loop(Sys, controller);

%%
% We can run our system in closed loop, but this is not very interesting,
% because w is 0 anyway. Let's change that 
Sys.Wref = Sys.time*0.;
Sys.Wref(30:40) = 1; 
Sys.Wref(60:80) = -0.5; 

%%
% and the specs:
Sys.stl_list = {'alw (ev_[0,1.] alw_[0,0.5] ( abs(y1(t)-w1(t)) < 0.1))'};
controller = get_controller(Sys);

%%
% This time we will only plot input and outputs, i.e., disable the state
% plotting:
Sys.plot_x =[]; % default was Sys.plot_x = [1 2]
run_deterministic(Sys, controller);

%%
% More examples are given in the folder BluSTL/examples. In particular the 
% hvac_room case study demonstrate the adversarial scenario, as well as
% plot customization. The idea is to create a class derived from STLC_lti
% and specialize the update_plot method. 

