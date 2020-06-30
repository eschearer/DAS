%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HoldPositionAngleFinal.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MATLAB
clear
close all

tic;

%% Set up the simulation

%%% To improve stability, tighten the joint limits. A stiffness is applied
%%% when outside these limits
load model_struct
% Set scapula limits
model.dofs{1}.range=[-0.419 -0.419];
model.dofs{2}.range=[0.0874 0.0874];
model.dofs{3}.range=[12 12]*pi/180;
% Set AC limits
model.dofs{4}.range=[34 34]*pi/180;
model.dofs{5}.range=[-20 -20]*pi/180;
model.dofs{6}.range=[-16 -16]*pi/180;

% Initialize the model
das3('Initialize',model);
disp('Done.');

% model variables
ndof=11;
nmus=138;
nstates=2*ndof+2*nmus; % [Joint pos; joint vel; muscle LCE;muscle activation]

% indices to the state variables within the state vector
iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);
iAct = max(iLce) + (1:nmus);


%%% Set initial state
% The feasible states were determined by dragging the hand to the position using a
% robot and recording the final position and saved in 'feasiblepoints.mat'
pos = 8078; % select desired position index from feasible points
% Set initial state
load('feasiblepoints.mat')
x=zeros(nstates,1);
x(1:11)=stateFeasible(:,pos);

%%% We needed to pick a reasonable starting LCE for the muscles
% Initialize LCE to be right at slack length so force occurs immediately
% when the muscle activates.
LCEopt=das3('LCEopt');
muscle_tendon_lengths = das3('Musclelengths', x);
slack_lengths = das3('SEEslack');
Lce = muscle_tendon_lengths - slack_lengths;
x(iLce)=Lce; 

% Starting hand position
Phand_start=wristFeasible(pos,:)';

%%% Goal position is the position to move to. You can make this a
%%% trajectory and step through the goal in the simulation loop to move the
%%% hand. This code holds the wrist in a static position so the goal is the
%%% start.
% Set goal position
HandGoal=Phand_start;

%% Simulate the motion
complete=0;
while complete==0
    % Initialize derivatives and muscle excitations
    xdot = zeros(nstates,1);
    step_u = zeros(nmus,1);
    
    % Set simulation parameters
    time = 0;
    tend = 3;
    tstep = .003; % sometimes have to lessen this for stability
    nsteps = round((tend-time)/tstep)+1; % +1 allows for saving of initial state
    
    
    % initialize variables for saving
    Fhand=zeros(nsteps,3);
    usave=zeros(nsteps,nmus);
    xsave=zeros(nsteps,length(x));
    HandLocation=zeros(nsteps,3);
    JointAngles=zeros(nsteps,5);
    FBForce=zeros(nsteps,3);
    FBTorque=zeros(nsteps,4);
    TotalTorque=zeros(nsteps,4);
    Fkstep=zeros(nsteps,3); % proportional force called for
    Fdstep=zeros(nsteps,3); % derivative torque called for
    Fistep=zeros(nsteps,3); % integral force called for
    actStep=zeros(nsteps,9);
    
    
    
    % Initialize parameters
    u=zeros(nmus,1);
    handF=[0;0;0];
    
    % Save initial points
    i=1;
    Fhand(i,:)=handF';
    usave(i,:)=u';
    xsave(i,:)=x';
    HandLocation(i,:)=Phand_start';
    JointAngles(i,:)=qFeasible(pos,:);
    
    staticTorque=torqueFeasible(pos,:);
    
    M=Mfeasible(:,:,pos);
    
    openLoopAct=activationFeasible(pos,:);
    
    fprintf('\nSimulating...        ')
    sumErr=0;
    sumErrTau=0;
    alpha0=zeros(9,1);    
    
    %% Simulation loop
    % Run simulation
    for i=2:nsteps
        lastwarn('');
        a=0;
        % Display Progress
        if round(1000*i/nsteps)/10 < 10
            fprintf('\b\b\b\b\b%3.1f%%\n', round(1000*i/nsteps)/10)
        elseif round(1000*i/nsteps)/10 < 100
            fprintf('\b\b\b\b\b\b%3.1f%%\n', round(1000*i/nsteps)/10)
        elseif round(1000*i/nsteps)/10 < 1000
            fprintf('\b\b\b\b\b\b\b%3.1f%%\n', round(1000*i/nsteps)/10)
        end
        
        Moments=zeros(5,1);
        exF= [0;0]; % arm support of only a vertical force
        
        % Arm Support applied as stiffness and damping forces at the wrist
        [dPhand_dx, Phand] = pos_jacobian(x,model);
        supportEq=[0;0.3;-0.15];
        K=diag([0 30 30]);
        B=diag([120 120 120]);
        Vhand=dPhand_dx*x(12:22);
        handF=-K*(Phand-supportEq)-B*Vhand;
        
        %%% Apply feedback controller to the arm
        % Feedback parameters for the controller
        Kp=250;
        Kd=0;
        Ki=80;
        
        % Calculate feedback force to move wrist to desired position
        sumErr=sumErr+( HandGoal-Phand)*tstep;
        Fk=-Kp*(Phand-HandGoal);
        Fd=-Kd*Vhand;
        Fi=Ki*sumErr;
        Fdes=Fk+Fd+Fi;  
        
        % Convert feedback force to torque
        [~, ~, ~, ~, ~, ~, qTH] = das3('Dynamics',x, zeros(size(x)), zeros(138,1));
        pose=[qTH;x(10:11)];
        
        J=computeJacobianDAS_5angles(pose);
        J=J(:,1:4);
        FBtau=J'*Fdes;     
        
        tauDes=FBtau+staticTorque';
        
        % Holds the arm in place for 166 time steps
        % Goal is to get the simulation in a reasonable start state
        holdTime=166;
        if i<holdTime
            Khold=eye(3)*2000;
            Bhold=eye(3)*200;
            holdF=-Khold*(Phand-Phand_start)-Bhold*Vhand;
            handF=handF+holdF;
            K=K+Khold;
            B=B+Bhold;
        end
        
        %%% Compute and apply muscle activations. The activations can only
        %%% update at the stimulation frequency.
        if mod(i,25)==0 % only update at stimulation frequency of 13 Hz
            alpha0=computeActivations(M,tauDes,alpha0);
            for j=1:9
                mus=whichMuscles(j);
                u(mus)=alpha0(j)*1;
            end
        end
        %%% Advance simulation by a step
        [x, xdot, step_u] = das3step(x, u, tstep, xdot, step_u, Moments, exF, handF, K, B);
        
        
        % Save variables
        Fhand(i,:)=handF;
        usave(i,:)=u;
        xsave(i,:)=x;
        HandLocation(i,:)=Phand;
        JointAngles(i,:)=pose;
        FBForce(i,:)=Fdes;
        FBTorque(i,:)=FBtau;
        TotalTorque(i,:)=tauDes;
        Fkstep(i,:)=Fk;
        Fdstep(i,:)=Fd;
        Fistep(i,:)=Fi;
        actStep(i,:)=alpha0;
        
        
        %%% Catch warnings which means the answer is blowing up
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            a=1;
            break;
        end
    end
    %%% If a warning occurs, try a new LCE starting point and see if it
    %%% works
    if a==1
        complete=0;
        x=zeros(nstates,1);
        x(1:11)=stateFeasible(:,pos);
        x(iLce)=-3+6*rand(nmus,1); % randomly select initial Lce
    else
        complete=1; % simulation completed with no warnings
    end
end

%% Results and plots

% Final error
endErr=HandGoal'-HandLocation(end,:)
err=pdist([HandGoal';HandLocation(end,:)])


tout=tstep*(0:nsteps-1);
figure(1)
plot(tout,HandLocation);
title('Hand Position')

figure(2)
%plot(tout,Fhand)
plot(tout,FBTorque)
title('FB Torque')


figure(3)
plot(tout,actStep)
title('Activations')

figure(4)
plot(tout,JointAngles)
title('Joint Angles')

make_osimm('simulation', xsave(:,1:ndof), tout);

toc














