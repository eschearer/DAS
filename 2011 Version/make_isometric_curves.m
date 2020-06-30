function exitflag = make_isometric_curves
%Finds the maximum moment for various postures of the arm model,
% ensuring glenohumeral and scapular stability

% directory with output files
mydir = 'opt_results/';

params.alg = 'interior-point';
params.maxfulevals = 1000;
params.TolX = 1e-6;
params.TolFun = 1e-2;
params.TolCon = 1e-2;
% solver: either fmincon or ipopt
params.solver = 'fmincon';

% model sizes
ndof = 11;
nmus = 138;
nstates = 2*ndof + 2*nmus;

% indices to the state variables within the state vector
iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);
iAct = max(iLce) + (1:nmus);

aphi = 0.6728;      % 38.55 degrees
atheta = 0.7744;    % 44.37 degrees

aphitan=tand(38.55);
athetatan=tand(44.37);
colours = [0 0 0;247 147 30]/255;
set(0,'DefaultAxesFontSize',10);

% which dofs?
% (1)  elbow extension
% (2)  elbow flexion
% (3)  forearm supination
% (4)  forearm pronation
% (5)  shoulder extension
% (6)  shoulder flexion
% (7)  shoulder adduction
% (8)  shoulder abduction
% (9)  shoulder external rotation
% (10) shoulder internal rotation

all_dofs = [0 0 0 0 0 0 0 0 0 0];
%all_dofs = [1 1 1 1 1 1 1 0 1 1];

% which curves to plot?
%plot_dofs = [0 0 0 0 0 0 0 0 0 0];
plot_dofs = [0 0 0 0 0 0 0 0 1 1];
abstract = 0;

%% elbow flexion moment
elbow_flexion = [5 25 50 75 100 125]*pi/180;
momdof = 4;
momsign = 1;
fileindex=1;

% [0 pi/2 0 elbow_flexion(i) pi/2]: Buchanan, Amis
% [0 0 0 elbow_flexion(i) 0]: Garner, Winters

% maximise extension
for i=1:6
    arm_pos = [0 pi/2 0 elbow_flexion(i) pi/2]+ 1e-3;
    if all_dofs(1), run_opt; end
    fileindex=fileindex+1;
end
% maximise flexion
momsign = -1;
for i=1:6
    arm_pos = [0 pi/2 0 elbow_flexion(i) pi/2]+ 1e-3;
    if all_dofs(2) && i==5, run_opt; end
    fileindex=fileindex+1;
end

%% pronation/supination moment
pro_sup = [5 25 50 75 105 125]*pi/180;
momdof = 5;
momsign = 1;

%[0 0 0 pi/2 pro_sup(i)]: Garner, Salter, Winters
%[0 pi/2 0 pi/2 pro_sup(i)]: Winters

% maximise supination
for i=1:6
    arm_pos = [0 pi/18 0 pi/2 pro_sup(i)]+ 1e-3;
    if all_dofs(3) && i>3, run_opt; end
    fileindex=fileindex+1;
end
% maximise pronation
momsign = -1;
for i=1:6
    arm_pos = [0 pi/18 0 pi/2 pro_sup(i)]+ 1e-3;
    if all_dofs(4), run_opt; end
    fileindex=fileindex+1;
end

%% shoulder flexion/extension moment
shoulder_flexion = [5 25 40 55 70 80];
momdof = 2;
momsign = 1;

%[pi/2 shoulder_flexion(i) -pi/2 pi/3 pi/2]: Garner, Otis?
%[pi/2 shoulder_flexion(i) -pi/2 pi/2 pi/2]: Winters

% maximise extension
for i=1:6
    arm_pos = [85 shoulder_flexion(i) -85 60 90]*pi/180 + 1e-3;
    if all_dofs(5), run_opt; end
    fileindex=fileindex+1;
end
% maximise flexion
momsign = -1;
for i=1:6
    arm_pos = [85 shoulder_flexion(i) -85 60 90]*pi/180 + 1e-3;
    if all_dofs(6), run_opt; end
    fileindex=fileindex+1;
end

%% shoulder abduction moment
shoulder_abduction = [5 20 40 55 75 90]*pi/180;
momdof = 2;
momsign = 1;

%[pi/4 shoulder_abduction(i) 0 pi/3 pi/2]: Garner, Otis?

% maximise adduction
for i=1:6
    arm_pos = [pi/4 shoulder_abduction(i) 0 pi/3 pi/2]+ 1e-3;
    if all_dofs(7) && i==3, run_opt; end
    fileindex=fileindex+1;
end
% maximise abduction
momsign = -1;
for i=1:6
    arm_pos = [pi/4 shoulder_abduction(i) 0 pi/3 pi/2]+ 1e-3;
    if all_dofs(8), run_opt; end
    fileindex=fileindex+1;
end

%% shoulder rotation moment
shoulder_rotation = [-90 -65 -45 -30 -15 0]*pi/180;
momdof = 3;
momsign = 1;

%[pi/4 pi/3 shoulder_rotation(i) pi/3 pi/2]: Garner, Otis?

% maximise external rotation
for i=1:6
    arm_pos = [pi/4 pi/3 shoulder_rotation(i) pi/3 pi/2]+ 1e-3;
    if all_dofs(9) && i==1, run_opt; end
    fileindex=fileindex+1;
end
% maximise internal rotation
momsign = -1;
for i=1:6
    arm_pos = [pi/4 pi/3 shoulder_rotation(i) pi/3 pi/2]+ 1e-3;
    if all_dofs(10), run_opt; end
    fileindex=fileindex+1;
end

plot_exp_model;
%plot_GH;

% figure(1);
% export_fig elbow_fl_mom -pdf
% figure(2);
% export_fig forearm_pro_mom -pdf
% figure(3);
% export_fig sh_fl_mom -pdf
% figure(4);
% export_fig sh_abd_mom -pdf
% figure(5);
% export_fig sh_rot_mom -pdf

    function run_opt
        % Finds the maximum moment for a given posture of the arm model,
        % ensuring glenohumeral and scapular stability
        %
        % Calls "isometric_strength_sim" repeatedly if no solution is found, and
        % eventually saves the outputs from each call in opt_results\<optim_#>
        
 
        prev_res = load([mydir 'optim_51']);
        %prev_res = load([mydir 'optim_' num2str(fileindex) 'b']);
        state = [prev_res.osim_angles(:,end); prev_res.all_vel(:,end); prev_res.all_Lce(:,end); prev_res.all_act(:,end)];
        moments = prev_res.all_moments(:,end);
        x_prev = [state; moments];
        noprev=0;
        clear prev_res;
        
%        noprev=1;

%         exitflag=-10;
% %         prev_res = load([mydir 'optim_' num2str(fileindex)]);
% %         exitflag = prev_res.exitflag;
% %         clear prev_res;
%         while ~exitflag || exitflag==-10
%             if exitflag~=-10
%                 % this has already run once, but has not yet reached a
%                 % solution --> continue the search
%                 prev_res = load([mydir 'optim_' num2str(fileindex)]);
%                 state = [prev_res.osim_angles(:,end); prev_res.all_vel(:,end); prev_res.all_Lce(:,end); prev_res.all_act(:,end)];
%                 moments = prev_res.all_moments(:,end);
%                 x_prev = [state; moments];
%                 noprev=0;
%                 clear prev_res;
%             else
%                 % this is running for the first time
%                 noprev=1;
%                 pos6 = mod(fileindex,6);
%                 if pos6==0, pos6=6; end
%                 for filei=1:pos6-1
%                     prev_res = load([mydir 'optim_' num2str(fileindex-filei)]);
%                     if prev_res.exitflag~=0
%                         state = [prev_res.osim_angles(:,end); prev_res.all_vel(:,end); prev_res.all_Lce(:,end); prev_res.all_act(:,end)];
%                         moments = prev_res.all_moments(:,end);
%                         x_prev = [state; moments];
%                         noprev=0;
%                         clear prev_res;
%                         break;
%                     end
%                     clear prev_res;
%                 end
%             end
                        
            disp(['File optim_' num2str(fileindex) '...']);
            if noprev
                [state,moments,~,exitflag,output] = isometric_strength_sim(momdof,momsign,arm_pos,params);
            else
                [state,moments,~,exitflag,output] = isometric_strength_sim(momdof,momsign,arm_pos,params,x_prev);
            end
            
            x_prev = [state; moments];
            % save angles
            osim_angles = state(iq);
            % save angular velocities
            all_vel = state(iqdot);
            % save Lce
            all_Lce = state(iLce);
            % save activations
            all_act = state(iAct);
            % save forces
            all_forces = das3mex('Muscleforces', state);
            [~,opt_pos1] = simm2dsem(state(1:3),state(4:6),state(7:9));
            % save thoraco-humeral angles
            opt_pos = opt_pos1;
            % save moments
            all_moments = moments;
            % save scapular constraint
            all_Fscap = das3mex('Scapulacontact', state);
            FGH = das3mex('GHforce', x_prev(1:nstates));
            all_FGH = FGH;
            Rgt = glenoid_scap;
            Fgh0 = Rgt*FGH;  % take glenoid orientation into account
            if norm(Fgh0), Fgh0 = Fgh0/norm(Fgh0); end
            % decompose into polar angles
            thetar = asin(-Fgh0(2));
            if ~(sqrt(Fgh0(1)^2+Fgh0(3)^2)), phir = 0.0;
            else phir=asin(Fgh0(3)/sqrt(Fgh0(1)^2+Fgh0(3)^2));
            end
            FGHcontact = (thetar/atheta)^2 + (phir/aphi)^2 - 1; % <=0
            % save glenohumeral constraint
            all_FGHcontact = FGHcontact;
            
            save([mydir 'optim_' num2str(fileindex)],'osim_angles','all_vel','all_Lce',...
                'all_moments','all_act','all_forces','all_FGHcontact','all_FGH','all_Fscap',...
                'params','arm_pos','momdof','momsign','opt_pos','exitflag','output');
            make_osimm([mydir 'opt_angles_' num2str(fileindex)],osim_angles);
%        end
    end

    function plot_exp_model
        % Compare to experimental data from Holzbaur et al.
        
        leg_clav = {'protraction';'elevation';'axial rotation'};
        leg_scap = {'protraction';'lateral rotation';'spinal tilt'};
        colours = [0 0 0;247 147 30]/255;

         set(0, 'DefaultAxesFontSize', 11, ...
         'DefaultAXesColor', 'none',...
         'DefaultFigureColor', 'none')

        % elbow extension
        if plot_dofs(1)
            Buch = [20 30 40 50 60 70 80 90 100 110 120 130;...
                -17 -22 -30 -35 -40 -45 -48 -50 -52 -50 -45 -40;...
                  7   7  12  15  15  15  13  12  11   8   5   3;...
                  8   9  12  15  15  15  15  13  10   8   5   3];
            
            %             Amis = [0 25 50 75 100 120;...
            %                 -32 -40 -48 -52 -50 -45];
            
            Amis = [0 30 60 90 120 145;...
                -35 -41 -51 -52 -46 -42];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index)]);
                extmom(index) = mdata.all_moments(4);
                extang(index) = mdata.osim_angles(10)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
            end
            figure(1); %subplot(3,1,1); 
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            h1 = errorbar(Buch(1,:),Buch(2,:),Buch(4,:),Buch(3,:),'-','Color',colours(1,:));  hold on;
            h2 = plot(Amis(1,:),Amis(2,:),':','Color',colours(1,:));
            h3 = plot(extang(flags>0),-extmom(flags>0),'-','Color',colours(1,:),'Linewidth',1.5);
            text(10,-70,'Extension','Color',colours(1,:),'FontSize',10);
            
            if ~plot_dofs(2)
                lh = legend('Buchanan','Amis','Model');
                set(lh,'Box','on','Location','East');
                ylabel('Elbow flexion moment (Nm)');
            end
            
            subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(extang(flags>0),clav_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(extang(flags>0),clav_angles(2,flags>0),':','Color',colours(1,:));
            plot(extang(flags>0),clav_angles(3,flags>0),'--','Color',colours(1,:));
            xlim([0 150]);
            ylabel('Clavicular angles (degrees)');
            legend('protraction','elevation','axial rotation');
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(extang(flags>0),scap_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(extang(flags>0),scap_angles(2,flags>0),':','Color',colours(1,:));
            plot(extang(flags>0),scap_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Scapular angles (degrees)');
            legend('protraction','lateral rotation','spinal tilt');
            xlim([0 150]);
            
            if ~plot_dofs(2)
                xlabel('Elbow flexion angle (degrees)');
            end
            
       end
        
        % elbow flexion
        if plot_dofs(2)
            Buch = [20 30 40 50 60 70 80 90 100 110 120 130;...
                    40 44 49 54 60 65 72 75  76  74  70  58;...
                    15 11  9  9  9  9  9  9   9  11  11  15;...
                    15 11  9  9  9  9  9  9   9  11  11  15];
            
            %             Amis = [0 25 50 75 100 120;...
            %                 55 65 73 73 62 48];
            
            Amis = [0 30 60 90 120 145;...
                53 68 76 70 51 35];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+6)]);
                flexmom(index) = mdata.all_moments(4);
                flexang(index) = mdata.osim_angles(10)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
            end
            
            figure(1); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            errorbar(Buch(1,:),Buch(2,:),Buch(3,:),Buch(4,:),'-','Color',colours(2,:)); hold on;
            plot(Amis(1,:),Amis(2,:),':','Color',colours(2,:));
            plot(flexang(flags>0),-flexmom(flags>0),'-','Color',colours(2,:),'Linewidth',1.5);
            text(10,75,'Flexion','Color',colours(2,:),'FontSize',10);
            
            lh = legend([h1,h2,h3],'Buchanan','Amis','Model');
            set(lh,'EdgeColor','w','Location','East');
            ylabel('Nm');
            title('A: Elbow flexion/extension moment','FontWeight','bold','Position',[5,102],'HorizontalAlignment','left');
            xlim([0 135]); 
 
            subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(extang(flags>0),clav_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(extang(flags>0),clav_angles(2,flags>0),':','Color',colours(2,:));
            plot(extang(flags>0),clav_angles(3,flags>0),'--','Color',colours(2,:));
            ylabel('degrees');
            title('B: Clavicular angles','FontWeight','bold','Position',[5,102],'HorizontalAlignment','left');
            legend(leg_clav,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlim([0 135]); ylim([-120 100]);
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(extang(flags>0),scap_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(extang(flags>0),scap_angles(2,flags>0),':','Color',colours(2,:)); 
            plot(extang(flags>0),scap_angles(3,flags>0),'--','Color',colours(2,:)); 
            title('C: Scapular angles','FontWeight','bold','Position',[5,72],'HorizontalAlignment','left');
            ylabel('degrees');
            legend(leg_scap,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlabel('Elbow flexion angle (degrees)');            
            xlim([0 135]); ylim([-40 70]);
        end
        
        % forearm supination
        if plot_dofs(3)
            Garner = [-90 -75 -60 -40 -30 -25 -18 0 25;...
                -2 -8 -9 -8.5 -12 -13 -14 -15 -13];
            Winters = [-90 -60 -30 0 30 60;...
                -3 -5.5 -7.2 -9.1 -10.2 -10.9];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+2*6)]);
                supmom(index) = mdata.all_moments(5);
                supang(index) = mdata.osim_angles(11)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
            end
            
            figure(2); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            h1 = plot(Garner(1,:),Garner(2,:),'--','Color',colours(1,:)); hold on;
            h2 = plot(Winters(1,:),Winters(2,:),'-','Color',colours(1,:));
            h3 = plot(supang(flags>0)-90,-supmom(flags>0),'-','Color',colours(1,:),'Linewidth',1.5);
            text(-70,-11,'Supination','Color',colours(1,:),'FontSize',10);
            xlim([-90 70]);
            
            if ~plot_dofs(4)
                ylabel('Pronation moment (Nm)');
                lh=legend('Garner','Model');
                set(lh,'Box','on','Location','East');
            end
            
            subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(supang(flags>0)-90,clav_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(supang(flags>0)-90,clav_angles(2,flags>0),':','Color',colours(1,:));
            plot(supang(flags>0)-90,clav_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Clavicular angles (degrees)');
            legend('protraction','elevation','axial rotation');
            xlim([-90 70]);
           
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(supang(flags>0)-90,scap_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(supang(flags>0)-90,scap_angles(2,flags>0),':','Color',colours(1,:));
            plot(supang(flags>0)-90,scap_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Scapular angles (degrees)');
            legend('protraction','lateral rotation','spinal tilt');
            xlim([-90 70]);
            
            if ~plot_dofs(4)
                xlabel('Pronation angle (degrees)');
            end
        end
        
        % forearm pronation
        if plot_dofs(4)
            Garner = [-90 -75 -60 -40 -30 -25 -18 0 25;...
                14 9.5 10 10.5 11 10 9 9.5 7];
            Winters = [-90 -60 -30 0 30 60;...
                11.5 10.8 10 8 6 3];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+3*6)]);
                promom(index) = mdata.all_moments(5);
                proang(index) = mdata.osim_angles(11)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
            end
            
            figure(2); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            plot(Garner(1,:),Garner(2,:),'--','Color',colours(2,:));  hold on;
            plot(Winters(1,:),Winters(2,:),'-','Color',colours(2,:));
            plot(proang(flags>0)-90,-promom(flags>0),'-','Color',colours(2,:),'Linewidth',1.5);
            text(-70,13,'Pronation','Color',colours(2,:),'FontSize',10);
            xlim([-90 40]); ylim([-16 15]);
            
            ylabel('Nm');
            title('A: Pronation/supination moment','FontWeight','bold','Position',[-85,16],'HorizontalAlignment','left');
            lh=legend([h1,h2,h3],'Garner','Winters','Model');
            set(lh,'EdgeColor','w','Location','East');

            subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(proang(flags>0)-90,clav_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(proang(flags>0)-90,clav_angles(2,flags>0),':','Color',colours(2,:));
            plot(proang(flags>0)-90,clav_angles(3,flags>0),'--','Color',colours(2,:));
            ylabel('degrees');
            title('B: Clavicular angles','FontWeight','bold','Position',[-85,37],'HorizontalAlignment','left');
            legend(leg_clav,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlim([-90 40]); ylim([-70 35]); 
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(proang(flags>0)-90,scap_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(proang(flags>0)-90,scap_angles(2,flags>0),':','Color',colours(2,:)); 
            plot(proang(flags>0)-90,scap_angles(3,flags>0),'--','Color',colours(2,:)); 
            title('C: Scapular angles','FontWeight','bold','Position',[-85,13],'HorizontalAlignment','left');
            ylabel('degrees');
            legend(leg_scap,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlabel('Pronation angle (degrees)');
            xlim([-90 40]); ylim([-15 12]);       
        end
        
        % shoulder extension
        if plot_dofs(5)
            Garner = [0 20 40 60 80 90;...
                -75 -85 -87 -90 -88 -87];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+4*6)]);
                extmom(index) = mdata.all_moments(2);
                extang(index) = mdata.opt_pos(2)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
            end
            
            figure(3); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            h1 = plot(Garner(1,:),Garner(2,:),'--','Color',colours(1,:)); hold on;
            h4 = plot(extang(flags>0),-extmom(flags>0),'-','Color',colours(1,:),'Linewidth',1.5);
            text(60,-105,'Extension','Color',colours(1,:),'FontSize',10);
            
            if ~plot_dofs(6)
                ylabel('Shoulder flexion moment (Nm)');
                lh=legend('Garner','Model');
                set(lh,'Box','on','Location','West');
            end
            
            subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(extang(flags>0),clav_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(extang(flags>0),clav_angles(2,flags>0),':','Color',colours(1,:));
            plot(extang(flags>0),clav_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Clavicular angles (degrees)');
            legend('protraction','elevation','axial rotation');
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(extang(flags>0),scap_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(extang(flags>0),scap_angles(2,flags>0),':','Color',colours(1,:));
            plot(extang(flags>0),scap_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Scapular angles (degrees)');
            legend('protraction','lateral rotation','spinal tilt');
            
            if ~plot_dofs(6)
                xlabel('Shoulder flexion angle (degrees)');
            end
        end
        
        % shoulder flexion
        if plot_dofs(6)
            Garner = [0 20 40 60 80 90;...
                115 100 80 85 75 79];
            
            Otis = [0 45 90;...
                93.8 80.6 75.2;...
                4.4 3.1 2.5;...
                4.4 3.1 2.5];
            
            Winters = [0 15 30 45 60 75 90;...
                      92 85 75 70 68 62 60;...
                       4  4  4  9  8  7  6;...
                       4  4  4  9  8  7  6];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+5*6)]);
                flexmom(index) = mdata.all_moments(2);
                flexang(index) = mdata.opt_pos(2)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
           end
            
            figure(3); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            plot(Garner(1,:),Garner(2,:),'--','Color',colours(2,:));  hold on;
            h2=errorbar(Otis(1,:),Otis(2,:),Otis(3,:),Otis(4,:),':','Color',colours(2,:));
            h3=errorbar(Winters(1,:),Winters(2,:),Winters(3,:),Winters(4,:),'-','Color',colours(2,:));
            plot(flexang(flags>0),-flexmom(flags>0),'-','Color',colours(2,:),'Linewidth',1.5);
            text(60,100,'Flexion','Color',colours(2,:),'FontSize',10);
            xlim([-0.5 90.5]);
            ylim([-130 120]);
            
            ylabel('Nm');
            title('A: Shoulder flexion/extension moment','FontWeight','bold','Position',[4.5,122],'HorizontalAlignment','left');
            lh=legend([h1,h2,h3,h4],'Garner','Otis','Winters','Model');
            set(lh,'EdgeColor','w','Location','East');

            subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(flexang(flags>0),clav_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(flexang(flags>0),clav_angles(2,flags>0),':','Color',colours(2,:));
            plot(flexang(flags>0),clav_angles(3,flags>0),'--','Color',colours(2,:));
            ylabel('degrees');
            title('B: Clavicular angles','FontWeight','bold','Position',[4.5,24],'HorizontalAlignment','left');
            legend(leg_clav,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            ylim([-65 22]);
            xlim([-0.5 90.5]);
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(flexang(flags>0),scap_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(flexang(flags>0),scap_angles(2,flags>0),':','Color',colours(2,:)); 
            plot(flexang(flags>0),scap_angles(3,flags>0),'--','Color',colours(2,:)); 
            title('C: Scapular angles','FontWeight','bold','Position',[4.5,16],'HorizontalAlignment','left');
            ylabel('degrees');
            legend(leg_scap,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlabel('Shoulder flexion angle (degrees)');
            ylim([-15 15]);
            xlim([-0.5 90.5]);
        end
        
        % shoulder adduction
        if plot_dofs(7)
            Garner = [10 20 40 60 65 80 90;...
                -50 -60 -75 -85 -90 -80 -75];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+6*6)]);
                addmom(index) = mdata.all_moments(2);
                addang(index) = mdata.opt_pos(2)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
           end
            
            figure(4); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            h1 = plot(Garner(1,:),Garner(2,:),'--','Color',colours(1,:)); hold on;
            h3 = plot(addang(flags>0),-addmom(flags>0),'-','Color',colours(1,:),'Linewidth',1.5);
            text(30,-115,'Adduction','Color',colours(1,:),'FontSize',10);

            % for abstract
            if abstract
                figure(6); %subplot(2,1,1);
                subplot('Position',[0.1 0.6 0.85 0.35]);
                plot(Garner(1,:),Garner(2,:),'--','Color',colours(1,:)); hold on;
                plot(addang(flags>0),-addmom(flags>0),'-','Color',colours(1,:),'Linewidth',1.5);
                text(30,-115,'Adduction','Color',colours(1,:),'FontSize',10);
            end
        
            if ~plot_dofs(8)
                ylabel('Shoulder abduction moment (Nm)');
                lh=legend('Garner','Model');
                set(lh,'Box','on','Location','West');
            end
            
            figure(4); subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(addang(flags>0),clav_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(addang(flags>0),clav_angles(2,flags>0),':','Color',colours(1,:));
            plot(addang(flags>0),clav_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Clavicular angles (degrees)');
            legend('protraction','elevation','axial rotation');
            xlim([0 85]);
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(addang(flags>0),scap_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(addang(flags>0),scap_angles(2,flags>0),':','Color',colours(1,:));
            plot(addang(flags>0),scap_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Scapular angles (degrees)');
            legend('protraction','lateral rotation','spinal tilt');
            xlim([0 85]);
            
            if ~plot_dofs(8)
                xlabel('Shoulder abduction angle (degrees)');
            end
        end
        
        % shoulder abduction
        if plot_dofs(8)
            Garner = [10 20 40 60 65 80 90;...
                100 85 75 77 78 80 75];
            
            Otis = [0 45 90;...
                71.4 57.9 57.1;...
                3.3 1.8 2.3;...
                3.3 1.8 2.3];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+7*6)]);
                abdmom(index) = mdata.all_moments(2);
                abdang(index) = mdata.opt_pos(2)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
            end
            
            figure(4); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            plot(Garner(1,:),Garner(2,:),'--','Color',colours(2,:));  hold on;
            h2=plot(Otis(1,:),Otis(2,:),':','Color',colours(2,:));
            plot(abdang(flags>0),-abdmom(flags>0),'-','Color',colours(2,:),'Linewidth',1.5);
            text(30,95,'Abduction','Color',colours(2,:),'FontSize',10);
            ylim([-130 110]);
            xlim([-0.5 90.5]);
            
            ylabel('Nm');
            title('A: Shoulder abduction/adduction moment','FontWeight','bold','Position',[4.5,112],'HorizontalAlignment','left');
            lh=legend([h1,h2,h3],'Garner','Otis','Model');
            set(lh,'EdgeColor','w','Location','East');

            % for abstract
            if abstract
                figure(6); %subplot(2,1,1);
                subplot('Position',[0.1 0.6 0.85 0.35]);
                h11 = plot(Garner(1,:),Garner(2,:),'--','Color',colours(2,:));  hold on;
                h21 = plot(Otis(1,:),Otis(2,:),':','Color',colours(2,:));
                h31 = plot(abdang(flags>0),-abdmom(flags>0),'-','Color',colours(2,:),'Linewidth',1.5);
                text(30,95,'Abduction','Color',colours(2,:),'FontSize',10);
                ylim([-130 110]);
                xlim([-0.5 90.5]);

                ylabel('Nm');
                xlabel('Shoulder abduction angle (degrees)');
                title('A: Shoulder abduction/adduction moment','FontWeight','bold','Position',[4.5,112],'HorizontalAlignment','left');
                lh=legend([h11,h21,h31],'Garner','Otis','Model');
                set(lh,'EdgeColor','w','Location','East');
            end
            
            figure(4); subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(abdang(flags>0),clav_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(abdang(flags>0),clav_angles(2,flags>0),':','Color',colours(2,:));
            plot(abdang(flags>0),clav_angles(3,flags>0),'--','Color',colours(2,:));
            ylabel('degrees');
            title('B: Clavicular angles','FontWeight','bold','Position',[4.5,87],'HorizontalAlignment','left');
            legend(leg_clav,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlim([-0.5 90.5]); ylim([-60 85]);
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(abdang(flags>0),scap_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(abdang(flags>0),scap_angles(2,flags>0),':','Color',colours(2,:)); 
            plot(abdang(flags>0),scap_angles(3,flags>0),'--','Color',colours(2,:)); 
            title('C: Scapular angles','FontWeight','bold','Position',[4.5,72],'HorizontalAlignment','left');
            ylabel('degrees');
            legend(leg_scap,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlabel('Shoulder abduction angle (degrees)');
            xlim([-0.5 90.5]); ylim([-30 70]);
        
        end
        
        % shoulder external rotation
        if plot_dofs(9)
            Garner = [-90 -80 -60 -40 -20 0;...
                -45 -44 -45 -48 -49 -47];
            
            Otis = [-45 0;...
                -35.7 -43.3;...
                1.6 1.5;...
                1.6 1.5];
            
            Engin = [-90 -60 -30 0;...
                -17 -28 -32 -30;...
                -22 -32 -42 -48;...
                -39 -55 -50 -32];
                
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+8*6)]);
                extmom(index) = mdata.all_moments(3);
                extang(index) = mdata.opt_pos(3)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
            end
            
            figure(5); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            h1 = plot(Garner(1,:),Garner(2,:),'--','Color',colours(1,:));  hold on;
            h2 = errorbar(Otis(1,:),Otis(2,:),Otis(3,:),Otis(4,:),':','Color',colours(1,:));
            h3 = errorbar(Engin(1,:),mean(Engin(2:4,:)),std(Engin(2:4,:)),std(Engin(2:4,:)),'-','Color',colours(1,:));
            h4 = plot(extang(flags>0),-extmom(flags>0),'-','Color',colours(1,:),'Linewidth',1.5);
            text(-80,-55,'External rotation','Color',colours(1,:),'FontSize',10);
            
            % for abstract
            if abstract
                figure(6); %subplot(2,1,2);
                subplot('Position',[0.1 0.1 0.85 0.35]);
                plot(Garner(1,:),Garner(2,:),'--','Color',colours(1,:));  hold on;
                errorbar(Otis(1,:),Otis(2,:),Otis(3,:),Otis(4,:),':','Color',colours(1,:));
                plot(extang(flags>0),-extmom(flags>0),'-','Color',colours(1,:),'Linewidth',1.5);
                text(-80,-55,'External rotation','Color',colours(1,:),'FontSize',10);
            end

            if ~plot_dofs(10)
                ylabel('Shoulder rotation moment (Nm)');
                lh=legend('Garner','Otis','Engin','Model');
                set(lh,'Box','on','Location','West');
            end

            figure(5); subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(extang(flags>0),clav_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(extang(flags>0),clav_angles(2,flags>0),':','Color',colours(1,:));
            plot(extang(flags>0),clav_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Clavicular angles (degrees)');
            legend('protraction','elevation','axial rotation');
            xlim([-91 1]);
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(extang(flags>0),scap_angles(1,flags>0),'-','Color',colours(1,:)); hold on;
            plot(extang(flags>0),scap_angles(2,flags>0),':','Color',colours(1,:));
            plot(extang(flags>0),scap_angles(3,flags>0),'--','Color',colours(1,:));
            ylabel('Scapular angles (degrees)');
            legend('protraction','lateral rotation','spinal tilt');
            xlim([-91 1]);
            
            if ~plot_dofs(10)
                xlabel('Shoulder rotation angle (degrees)');
            end
        end
        
        % shoulder internal rotation
        if plot_dofs(10)
            Garner = [-90 -80 -60 -40 -20 0;...
                65 70 76 77 75 71];
            
            Otis = [-45 0;...
                51.1 47.8;...
                2.0 1.7;...
                2.0 1.7];
            
            Engin = [-90 -60 -30 0;...
                45 39 39 40;...
                65 61 50 43;...
                77 67 61 59];
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+9*6)]);
                intmom(index) = mdata.all_moments(3);
                intang(index) = mdata.opt_pos(3)*180/pi;
                flags(index) = mdata.exitflag;
                clav_angles(:,index) = mdata.osim_angles(1:3)*180/pi;
                [hvec,~] = simm2dsem(mdata.osim_angles(1:3),mdata.osim_angles(4:6),mdata.osim_angles(7:9));
                scap_angles(:,index) = hvec*180/pi;
            end
            
            figure(5); %subplot(3,1,1);
            set(gcf,'Position',[102    15   568   667]);
            subplot('Position',[0.1 0.64 0.85 0.31]);
            plot(Garner(1,:),Garner(2,:),'--','Color',colours(2,:)); hold on;
            errorbar(Otis(1,:),Otis(2,:),Otis(3,:),Otis(4,:),':','Color',colours(2,:));
            errorbar(Engin(1,:),mean(Engin(2:4,:)),std(Engin(2:4,:)),std(Engin(2:4,:)),'-','Color',colours(2,:));
            plot(intang(flags>0),-intmom(flags>0),'-','Color',colours(2,:),'Linewidth',1.5);
            text(-80,80,'Internal rotation','Color',colours(2,:),'FontSize',10);
            xlim([-91 1]);
            ylim([-70 90]);
            
            ylabel('Nm');
            title('A: Shoulder internal/external rotation moment','FontWeight','bold','Position',[-85,92],'HorizontalAlignment','left');
            lh = legend([h1,h2,h3,h4],{'Garner';'Otis';'Engin';'Model'});
            set(lh,'EdgeColor','w','Location','East');

             % for abstract
            if abstract
                figure(6); %subplot(2,1,2);
                subplot('Position',[0.1 0.1 0.85 0.35]);
                h11 = plot(Garner(1,:),Garner(2,:),'--','Color',colours(2,:)); hold on;
                h21 = errorbar(Otis(1,:),Otis(2,:),Otis(3,:),Otis(4,:),':','Color',colours(2,:));
                h32 = plot(intang(flags>0),-intmom(flags>0),'-','Color',colours(2,:),'Linewidth',1.5);
                text(-80,80,'Internal rotation','Color',colours(2,:),'FontSize',10);
                xlim([-90 0]);
                ylim([-70 90]);

                ylabel('Nm');
                xlabel('Shoulder rotation angle (degrees)');
                title('B: Shoulder internal/external rotation moment','FontWeight','bold','Position',[-85,92],'HorizontalAlignment','left');
                lh = legend([h11,h21,h31],{'Garner';'Otis';'Model'});
                set(lh,'EdgeColor','w','Location','East');
            end
            
            figure(5); subplot('Position',[0.1 0.36 0.85 0.2]); %subplot(3,1,2);
            plot(intang(flags>0),clav_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(intang(flags>0),clav_angles(2,flags>0),':','Color',colours(2,:));
            plot(intang(flags>0),clav_angles(3,flags>0),'--','Color',colours(2,:));
            ylabel('degrees');
            title('B: Clavicular angles','FontWeight','bold','Position',[-85,102],'HorizontalAlignment','left');
            legend(leg_clav,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlim([-91 1]);
            
            subplot('Position',[0.1 0.08 0.85 0.2]); %subplot(3,1,3);
            plot(intang(flags>0),scap_angles(1,flags>0),'-','Color',colours(2,:)); hold on;
            plot(intang(flags>0),scap_angles(2,flags>0),':','Color',colours(2,:)); 
            plot(intang(flags>0),scap_angles(3,flags>0),'--','Color',colours(2,:)); 
            title('C: Scapular angles','FontWeight','bold','Position',[-85,72],'HorizontalAlignment','left');
            ylabel('degrees');
            legend(leg_scap,'Orientation','horizontal','Location','SouthEast','EdgeColor','w');
            xlabel('Shoulder rotation angle (degrees)');
            xlim([-91 1]); ylim([-50 70]);
            
        end
        
    end

    function plot_GH
        % elbow extension/flexion
        if 0
            %        if plot_dofs(1) && plot_dofs(2)
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index)]);
                [thetar_ex(index),phir_ex(index)] = findGHangles(mdata.all_FGH);
            end
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+6)]);
                [thetar_fl(index),phir_fl(index)] = findGHangles(mdata.all_FGH);
            end
            
            scatter_GH(phir_ex,thetar_ex,phir_fl,thetar_fl,colours(1,:),colours(2,:),1);
        end
        
        % forearm supination/pronation
        if 0
            %        if plot_dofs(3) && plot_dofs(4)
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+2*6)]);
                [thetar_sup(index),phir_sup(index)] = findGHangles(mdata.all_FGH);
            end
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+3*6)]);
                [thetar_pro(index),phir_pro(index)] = findGHangles(mdata.all_FGH);
            end
            
            scatter_GH(phir_sup,thetar_sup,phir_pro,thetar_pro,colours(1,:),colours(2,:),2);
        end
        
        % shoulder extension/flexcion
        if plot_dofs(5) && plot_dofs(6)
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+4*6)]);
                [thetar_ex(index),phir_ex(index)] = findGHangles(mdata.all_FGH);
            end
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+5*6)]);
                [thetar_fl(index),phir_fl(index)] = findGHangles(mdata.all_FGH);
            end
            
            scatter_GH(phir_ex,thetar_ex,phir_fl,thetar_fl,colours(1,:),colours(2,:),3);
        end
        
        % shoulder adduction/abduction
        if plot_dofs(7) && plot_dofs(8)
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+6*6)]);
                [thetar_add(index),phir_add(index)] = findGHangles(mdata.all_FGH);
            end
            
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+7*6)]);
                [thetar_abd(index),phir_abd(index)] = findGHangles(mdata.all_FGH);
            end
            
            scatter_GH(phir_add,thetar_add,phir_abd,thetar_abd,colours(1,:),colours(2,:),4);
        end
        
        % shoulder external/internal rotation
        if plot_dofs(9) && plot_dofs(10)
            for index=1:6
                mdata = load([mydir 'optim_' num2str(index+8*6)]);
                [thetar_ex(index),phir_ex(index)] = findGHangles(mdata.all_FGH);
            end
            
            index=1;
            for rotindex=1:6
                mdata = load([mydir 'optim_' num2str(rotindex+9*6)]);
                [thetar_in(index),phir_in(index)] = findGHangles(mdata.all_FGH);
                index=index+1;
            end
            
            scatter_GH(phir_ex,thetar_ex,phir_in,thetar_in,colours(1,:),colours(2,:),5);
        end
    end

    function scatter_GH(phir1,thetar1,phir2,thetar2,colour1,colour2,fig_h)
        figure(fig_h);
        axes;
        Plellips(0,0,aphitan,athetatan);
        %ylim([-50 50]);
        axis equal
        axis on
        box off
        xlabel('Inferior','fontsize',12);
        ylabel('Posterior','fontsize',12);
        hold on
        scatter(-tan(phir1),-tan(thetar1),20,'b','filled');
        scatter(-tan(phir1(1)),-tan(thetar1(1)),20,'g','filled');
        plot(-tan(phir1),-tan(thetar1),'b');
        scatter(-tan(phir2),-tan(thetar2),20,'r','filled');
        scatter(-tan(phir2(1)),-tan(thetar2(1)),20,'g','filled');
        plot(-tan(phir2),-tan(thetar2),'r');
        %         plot(-tan(phir1),-tan(thetar1),'-+b');
        %         plot(-tan(phir2),-tan(thetar2),'-^r');
        set(gca,'xtick',[],'ytick',[]);
        
        inset_size=0.25;
        main_fig = findobj(fig_h,'Type','axes');
        ax=get(main_fig(3),'Position');
        set(main_fig(1),'Position', [ax(1)+ax(3)-inset_size -1.6*ax(2)+ax(4)-inset_size inset_size inset_size]);
    end

    function [thetar,phir] = findGHangles(FGH)
        Rgt = glenoid_scap;
        
        FGHsize = size(FGH);
        if FGHsize(2)==3, FGH = FGH'; end
        Fgh0 = Rgt*FGH;  % take glenoid orientation into account
        if norm(Fgh0), Fgh0 = Fgh0/norm(Fgh0); end
        % decompose into polar angles
        thetar = asin(-Fgh0(2));
        
        if ~(sqrt(Fgh0(1)^2+Fgh0(3)^2)), phir = 0.0;
        else phir=asin(Fgh0(3)/sqrt(Fgh0(1)^2+Fgh0(3)^2));
        end
    end
end
