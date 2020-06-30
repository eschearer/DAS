function error = get_momentarms(joints,muscles,osimfilename,mydir)
% builds .mat files that contain all moment arms of the DAS3 model
% inputs: structures "joints" and "muscles", created by "opensim_get_parameters.m"
% output: error = 0 if successful, -1 if not
%
% Dimitra Blana, February 2012
% based on code written for Wendy Murray's hand model
% 07/06/2012: Instead of opensim_get_momarm, use opensim_get_dLdq to
% calculate moment arms

import org.opensim.modeling.*
Mod = Model(osimfilename);

default_dof = 0; % default value for all dofs
allvectxts = '';
for ijnt=1:length(joints)
    jnt_values{ijnt} = default_dof;
    jnt_vals(ijnt) = default_dof;
    if isempty(allvectxts)
        allvectxts = ['jnt_values{1,' num2str(ijnt) '}'];
    else
        allvectxts = [allvectxts ',jnt_values{1,' num2str(ijnt) '}'];
    end
end

shoulder_angles = DAS3_workspace('DAS2data');

% save the muscle lengths when all dofs are zero
zeromusl = opensim_get_length(Mod, jnt_vals);
zeromusl = zeromusl(13:end);
save([mydir '\zerolength'],'zeromusl');

% for each muscle...
for imus = 1:length(muscles)
    mus = muscles{1,imus};
    
    % give all dofs the default value
    for ijnt=1:length(joints)
        jnt_values{ijnt} = default_dof;
    end
    
    % if the muscle only crosses the elbow, go through the range
    if mus.dofnums(1)>9
        for ijnt = 1:length(mus.dofnums)
            onedof = mus.dofnums(ijnt);
            lims = joints{1,onedof}.limits;
    %        values = lims(1):10:lims(2);  % every 10 degrees
            values = lims(1):(lims(2)-lims(1))/4:lims(2);  % 5 steps
            jnt_values{onedof} = values*pi/180;
        end

        % alljnts contain all the combinations of joint angles for the range of
        % motion of this muscle
        % rows: # combinations of angles
        % columns: # angles
        try
        eval(['alljnts = combvec(' allvectxts ')'';']);
        catch err
            disp(err);
            disp(['Skipping muscle ',mus.name, ' ...sorry!']);
            continue;     
        end
    elseif mus.dofnums(end)<10
        % if the muscle only crosses the shoulder, use "shoulder_angles"
        alljnts = shoulder_angles;
    else
        % if the muscle crosses both shoulder and elbow, use
        % shoulder_angles + range of motion
        clear elbow
        for ijnt = 10:11
            lims = joints{1,ijnt}.limits;
    %        values = lims(1):10:lims(2);  % every 10 degrees
            values = lims(1):(lims(2)-lims(1))/5:lims(2);  % 6 steps
            elbow{ijnt-9} = values*pi/180;
        end

        % alljnts contain all the combinations of joint angles for the range of
        % motion of this muscle
        % rows: # combinations of angles
        % columns: # angles
        try
        alljnts = combvec(shoulder_angles(:,1:9)',elbow{1},elbow{2})';
        catch err
            disp(err);
            disp(['Skipping muscle ',mus.name, ' ...sorry!']);
            continue;     
        end
    end
    
    % allmomarms contains the moment arms of this muscle for all
    % combinations of joint angles in alljnts
    % rows: # combinations of angles
    % columns: # of dofs this muscle crosses    
    allmomarms = zeros(size(alljnts,1),length(mus.dofnums));
    
    for ijnt = 1:length(mus.dofnums)
        all_dof_names{ijnt} = mus.dofnames{ijnt};
    end
    try
        % instead of moment arms directly from Opensim, calculate them using -dL/dq
        allmomarms = opensim_get_dLdq(Mod, alljnts, mus.name, all_dof_names);
        save([mydir '\momarms_' mus.name],'alljnts','allmomarms','jnt_values');
        disp([mus.name, ' moment arms saved.']);
        make_mot_file([mydir '\angles_' mus.name '.mot'],alljnts);
        disp(['Opensim motion file for muscle ', mus.name, ' created.']);
    catch err
        disp(err);
        error = -1;
        return;
    end

    clear mus alljnts allmomarms jnt_values all_dof_names
end

error = 0;