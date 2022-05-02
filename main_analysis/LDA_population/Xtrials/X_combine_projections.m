% load file
includedanimals = [1 2 3 4 5 6];
i0= 5;

proj_bs_go=cell(1,8);
proj_ls_go=cell(1,8);
proj_bs_nogo=cell(1,8);
proj_ls_nogo=cell(1,8);

for i = 1:8
    proj_bs_go{i} = [];
    proj_ls_go{i} = [];
    proj_bs_nogo{i} = [];
    proj_ls_nogo{i} = [];
    
    for animalnum =  includedanimals
        
        % trial per timepoint, averaged over reps
        alltrcrep_bs_go = [];
        alltrcrep_ls_go = [];
        alltrcrep_bs_nogo = [];
        alltrcrep_ls_nogo = [];
        % reps need to be averaged.
        for repi = 1: animalmodel{animalnum}.lmodel.nrep
            % 0*numtimepoints
            rep_bs_go = nan(0,size(animalmodel{animalnum}.lmodel.rep{1}.part{1}.go.alltraces{i},3));
            rep_ls_go = nan(0,size(animalmodel{animalnum}.lmodel.rep{1}.part{1}.go.alltraces{i},3));
            rep_bs_nogo = nan(0,size(animalmodel{animalnum}.lmodel.rep{1}.part{1}.nogo.alltraces{i},3));
            rep_ls_nogo = nan(0,size(animalmodel{animalnum}.lmodel.rep{1}.part{1}.nogo.alltraces{i},3));
            %%%%% 
            % numtrials in part 1 and 2 : both bs and ls
            % index is within each part
            trindbs_part1_go = find(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.Y{i} == 0);
            trindls_part1_go = find(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.Y{i} == 1);
            trindbs_part2_go = find(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.Y{i} == 0);
            trindls_part2_go = find(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.Y{i} == 1);
            
            trindbs_part1_nogo = find(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.Y{i} == 0);
            trindls_part1_nogo = find(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.Y{i} == 1);
            trindbs_part2_nogo = find(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.Y{i} == 0);
            trindls_part2_nogo = find(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.Y{i} == 1);
            
            % bs go
            % p1:vec p2:trials
            rep_bs_go = cattorep(rep_bs_go,trindbs_part1_go',i0,animalmodel{animalnum}.lmodel.rep{repi},2,1,i,'go');
            rep_bs_go = cattorep(rep_bs_go,trindbs_part2_go',i0,animalmodel{animalnum}.lmodel.rep{repi},1,2,i,'go');
            % ntrials*ntimepoints*nreps
            alltrcrep_bs_go = cat(3,alltrcrep_bs_go,rep_bs_go);
            
            % ls go
            % p1:vec p2:trials
            rep_ls_go = cattorep(rep_ls_go,trindls_part1_go',i0,animalmodel{animalnum}.lmodel.rep{repi},2,1,i,'go');
            rep_ls_go = cattorep(rep_ls_go,trindls_part2_go',i0,animalmodel{animalnum}.lmodel.rep{repi},1,2,i,'go');
            % ntrials*ntimepoints*nreps
            alltrcrep_ls_go = cat(3,alltrcrep_ls_go,rep_ls_go);
            
            % bs nogo
            % p1:vec p2:trials
            rep_bs_nogo = cattorep(rep_bs_nogo,trindbs_part1_nogo',i0,animalmodel{animalnum}.lmodel.rep{repi},2,1,i,'nogo');
            rep_bs_nogo = cattorep(rep_bs_nogo,trindbs_part2_nogo',i0,animalmodel{animalnum}.lmodel.rep{repi},1,2,i,'nogo');
            % ntrials*ntimepoints*nreps
            alltrcrep_bs_nogo = cat(3,alltrcrep_bs_nogo,rep_bs_nogo);
            
            % ls nogo
            % p1:vec p2:trials
            rep_ls_nogo = cattorep(rep_ls_nogo,trindls_part1_nogo',i0,animalmodel{animalnum}.lmodel.rep{repi},2,1,i,'nogo');
            rep_ls_nogo = cattorep(rep_ls_nogo,trindls_part2_nogo',i0,animalmodel{animalnum}.lmodel.rep{repi},1,2,i,'nogo');
            % ntrials*ntimepoints*nreps
            alltrcrep_ls_nogo = cat(3,alltrcrep_ls_nogo,rep_ls_nogo);
            
        end
        % per animal bias
        gobias = nanmean(nanmean([(squeeze(nanmean(alltrcrep_bs_go,3)));(squeeze(nanmean(alltrcrep_ls_go,3)))]));
        nogobias = nanmean(nanmean([(squeeze(nanmean(alltrcrep_bs_nogo,3)));(squeeze(nanmean(alltrcrep_ls_nogo,3)))]));
        % average over reps : ntrials*ntimepoints
        proj_bs_go{i} = cat(1,proj_bs_go{i},(squeeze(nanmean(alltrcrep_bs_go,3)))-gobias);
        proj_ls_go{i} = cat(1,proj_ls_go{i},(squeeze(nanmean(alltrcrep_ls_go,3)))-gobias);
        
        proj_bs_nogo{i} = cat(1,proj_bs_nogo{i},(squeeze(nanmean(alltrcrep_bs_nogo,3)))-nogobias);
        proj_ls_nogo{i} = cat(1,proj_ls_nogo{i},(squeeze(nanmean(alltrcrep_ls_nogo,3)))-nogobias);
    end
    
end

smoothspan = 5;
coef = 1
figure;
for i=1:8
    subplot(8,1,i);
    shadedErrorBar([],smooth(nanmean(proj_bs_go{i},1),smoothspan),...
        smooth(coef*nanstd(proj_bs_go{i})/sqrt(size(proj_bs_go{i},1)),smoothspan),'k',1);
    hold on;
    shadedErrorBar([],smooth(nanmean(proj_ls_go{i},1),smoothspan),...
        smooth(coef*nanstd(proj_ls_go{i})/sqrt(size(proj_ls_go{i},1)),smoothspan),'b',1)
    
end

figure;
for i=1:8
    subplot(8,1,i);
    shadedErrorBar([],smooth(nanmean(proj_bs_nogo{i},1),smoothspan),...
        smooth(coef*nanstd(proj_bs_nogo{i})/sqrt(size(proj_bs_nogo{i},1)),smoothspan),'k',1);
    hold on;
    shadedErrorBar([],smooth(nanmean(proj_ls_nogo{i},1),smoothspan),...
        smooth(coef*nanstd(proj_ls_nogo{i})/sqrt(size(proj_ls_nogo{i},1)),smoothspan),'b',1)
    
end
