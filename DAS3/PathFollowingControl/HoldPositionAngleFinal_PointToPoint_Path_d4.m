%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HoldPositionAngleFinal.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MATLAB
%clear
close all

tic;
switchTime=1.5; % Check swith modulus

% Initialize the model
load model_struct
model.dofs{1}.range=[-0.419 -0.419];
model.dofs{2}.range=[0.0874 0.0874];
model.dofs{3}.range=[12 12]*pi/180;
model.dofs{4}.range=[34 34]*pi/180;
model.dofs{5}.range=[-20 -20]*pi/180;
model.dofs{6}.range=[-16 -16]*pi/180;
das3('Initialize',model);
disp('Done.');

% model variables
ndof=11;
nmus=138;
nstates=2*ndof+2*nmus;

% indices to the state variables within the state vector
iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);
iAct = max(iLce) + (1:nmus);

d=4;
reachdist=30;



% Select target position
for trials=192:205
    for switchTime=[1.3]
        for D=[4]%,6,100]
            display(['D = ',num2str(d),' Trial: ',num2str(trials),' SwitchTime: ',num2str(switchTime)])
            filename=['Trial',num2str(trials),'Distance',num2str(100),'Time',num2str(floor(switchTime)),'.mat'];
            load(filename)
            paths
            d=D;
            if d==6 && switchTime == 1.3
                paths=FindPath(d,reachdist);
                paths=[paths paths(end)];
            elseif d==100
                pathlength=length(paths);
                pathStart=paths(1);
                pathEnd=paths(end);
                paths=pathStart;
                for i=1:pathlength-1
                    paths=[paths pathEnd];
                end
            else
                paths=FindPath(d,reachdist,paths(1),paths(end));
                paths=[paths paths(end)];
            end
            paths
            
            % Set initial state
            load('feasiblepoints.mat')
            pos=paths(1);
            x=zeros(nstates,1);
            x(1:11)=stateFeasible(:,pos);
            
            LCEopt=das3('LCEopt');
            muscle_tendon_lengths = das3('Musclelengths', x);    % only the first 11 elements of x (the joint angles) will be used
            slack_lengths = das3('SEEslack');
            Lce = muscle_tendon_lengths - slack_lengths;
            x(iLce)=Lce; % randomly select initial Lce
            
            Phand_start=wristFeasible(pos,:)';
            
            
            
            count=0;
            warnCount=0;
            while count==0
                % Set goal position
                pathIdx=2;
                GoalPos=paths(pathIdx);
                HandGoal=wristFeasible(paths(pathIdx),:)';
                % Initialize derivatives and muscle excitations
                xdot = zeros(nstates,1);
                step_u = zeros(nmus,1);
                
                % Set simulation parameters
                time = 0;
                tend = switchTime*(length(paths)-1);
                tstep = .003; % default was 0.003; decreased to make more stable numerically
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
                Fkstep=zeros(nsteps,3);
                Fdstep=zeros(nsteps,3);
                Fistep=zeros(nsteps,3);
                actStep=zeros(nsteps,9);
                GoalLocation=zeros(nsteps,3);
                
                
                
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
                GoalLocation(i,:)=HandGoal;
                
                staticTorque=torqueFeasible(GoalPos,:);
                
                M=Mfeasible(:,:,GoalPos);
                
                openLoopAct=activationFeasible(GoalPos,:);
                
                fprintf('\nSimulating...        ')
                sumErr=0;
                sumErrTau=0;
                alpha0=zeros(9,1);
                
                
                % % set first step activation
                % alpha0=computeActivations(M,staticTorque',openLoopAct');
                % for j=1:9
                %     mus=whichMuscles(j);
                %     u(mus)=alpha0(j)*1;
                % end
                % actStep(1,:)=alpha0;
                
                timer=0;
                % Run simulation
                for i=2:nsteps
                    timer=timer+tstep;
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
                    
                    if mod(i,floor(switchTime/tstep))==0 && i<nsteps-1.5*floor(switchTime/tstep) % makes sure last steps don't switch
                        pathIdx=pathIdx+1;
                        GoalPos=paths(pathIdx);
                        HandGoal=wristFeasible(paths(pathIdx),:)';
                        staticTorque=torqueFeasible(GoalPos,:);
                        M=Mfeasible(:,:,GoalPos);
                        openLoopAct=activationFeasible(GoalPos,:);
                        if d~=100
                            sumErr=0;
                        end
                    end
                    Moments=zeros(5,1);
                    exF= [0;0];
                    
                    % Arm Support
                    [dPhand_dx, Phand] = pos_jacobian(x,model);
                    supportEq=[0;0.3;-0.15];
                    K=diag([0 30 30]);
                    B=diag([120 120 120]);
                    Vhand=dPhand_dx*x(12:22);
                    handF=-K*(Phand-supportEq)-B*Vhand;
                    
                    holdTime=166;
                    if i<holdTime
                        Khold=eye(3)*2000;
                        Bhold=eye(3)*200;
                        holdF=-Khold*(Phand-Phand_start)-Bhold*Vhand;
                        handF=handF+holdF;
                        
                        
                        K=K+Khold;
                        B=B+Bhold;
                        
                    end
                    
                    % Feedback
                    Kp=250;
                    Kd=0;
                    Ki=80;
                    
                    sumErr=sumErr+(HandGoal-Phand)*tstep;
                    Fk=-Kp*(Phand-HandGoal);
                    Fd=-Kd*Vhand;
                    Fi=Ki*sumErr;
                    Fdes=Fk+Fd+Fi;
                    
                    [~, ~, ~, ~, ~, ~, qTH] = das3('Dynamics',x, zeros(size(x)), zeros(138,1));
                    pose=[qTH;x(10:11)];
                    
                    J=computeJacobianDAS_5angles(pose);
                    J=J(:,1:4);
                    FBtau=J'*Fdes;
                    
                    %     error=qFeasible(pos,:)'-pose;
                    %     sumErrTau=sumErrTau+error*tstep;
                    %     FBtau=Kp*error+Ki*sumErrTau;
                    %     FBtau=FBtau(1:4);
                    
                    
                    tauDes=FBtau+staticTorque';
                    %tauDes=[0;0;0;5];
                    
                    if mod(i,25)==0
                        alpha0=computeActivations(M,tauDes,alpha0);
                        for j=1:9
                            mus=whichMuscles(j);
                            u(mus)=alpha0(j)*1;
                        end
                    end
                    % Advance simulation by a step
                    [x, xdot, step_u] = das3step(x, u, tstep, xdot, step_u, Moments, exF, handF, K, B);
                    
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
                    GoalLocation(i,:)=HandGoal;
                    
                    [warnMsg, warnId] = lastwarn;
                    if ~isempty(warnMsg)
                        warnCount=warnCount+1;
                        a=1;
                        break;
                    end
                end
                if a==1 && warnCount<20
                    count=0;
                    x=zeros(nstates,1);
                    x(1:11)=stateFeasible(:,pos);
                    x(iLce)=-3+6*rand(nmus,1); % randomly select initial Lce
                else
                    count=1;
                end
            end
            
            endErr=HandGoal'-HandLocation(end,:)
            
            tout=tstep*(0:nsteps-1);
            
            mvar=mean(mean(MvarFeasible(:,:,GoalPos)));
            tvar=mean(mean(tvarFeasible(GoalPos,:)));
            
            filename=['Trial',num2str(trials),'Distance',num2str(d),'Time',num2str(floor(switchTime)),'.mat'];     
            
            figure(1)
            plot(tout,HandLocation);
            hold on
            plot(tout,GoalLocation);
            title('Hand Position')
            hold off
            %
            figure(2)
            %plot(tout,Fhand)
            plot(tout,FBTorque)
            title('FB Torque')
            
            drawnow
            
            
                    figure(3)
                    plot(tout,actStep)
                    title('Activations')
            
                    figure(4)
                    plot(tout,JointAngles)
                    title('Joint Angles')
            %
            
            make_osimm(filename, xsave(:,1:ndof), tout);
            
            toc
        end
    end
end












