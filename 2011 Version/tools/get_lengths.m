function error = get_lengths(joints,muscles,osimfilename,mydir)
% builds .mat files that contain all muscle-tendon lengths of the DAS3 model
% inputs: structures "joints" and "muscles", created by "opensim_get_parameters.m"
% output: error = 0 if successful, -1 if not
%
% Dimitra Blana, June 2012

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

%shoulder_angles = DAS3_workspace(-90:20:90,5:20:120,-55:20:70);
shoulder_angles = DAS3_workspace('DAS2data');

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
    
    % alllengths contains the length of this muscle for all
    % combinations of joint angles in alljnts
    % rows: # combinations of angles
    alllengths = zeros(size(alljnts,1),1);
    
    try
        allmuslength = opensim_get_length(Mod, alljnts);
        alllengths = allmuslength(:,imus+12);
        save([mydir '\lengths_' mus.name],'alljnts','alllengths','jnt_values');
        disp([mus.name, ' lengths saved.']);
    catch err
        disp(err);
        error = -1;
        return;
    end

    clear mus alljnts alllengths jnt_values
end

error = 0;