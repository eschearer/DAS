function change_momarms_bio(muscles,mydir)
% updates the bio file with new moment arm coefficients from a mat file

momarm_coef = load([mydir '/poly_results']);

biofile = fopen('../das3.bio','r');
newbiofile = fopen([mydir '/das3.bio'],'w');

for imus = 1:length(muscles)
    
    onemus = muscles{1,imus};
    num_musdofs = length(onemus.dofnums);
    musdofs = onemus.dofnums;
    
    mus = momarm_coef.mus_model{1,imus};
        
    line = fgetl(biofile);
    while (not(strncmp(line, ' parameters',11)))
        fprintf(newbiofile, '%s\r\n', line);
        line = fgetl(biofile);
    end
    num_params = sscanf(line, '%*s%i');
    for i=1:num_params
        line = fgetl(biofile);
    end
    
    fprintf(newbiofile, '%s%5i\n', ' parameters', mus.num_lparams);            
    for iparam = 1:mus.num_lparams
        for idof = 1:num_musdofs
            fprintf(newbiofile, '   %d ',mus.lparams(iparam,idof));            
        end        
        fprintf(newbiofile, '   %f\n',mus.lcoef(iparam));            
    end    
    
    fprintf('(%i) %s finished\n',imus,onemus.name);
    
end

% copy the rest of the file
while ~feof(biofile)
    line = fgetl(biofile);
    fprintf(newbiofile, '%s\r\n', line);
end

fclose(biofile);
fclose(newbiofile);
