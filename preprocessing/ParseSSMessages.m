%% behavior analysis

% 3: lick
% 4: running speed
% 7: beep


SSmessages = cell(size(expname));
correcttrial = cell(size(expname));
falsetrial = cell(size(expname));
misstrial = cell(size(expname));
alltrialcf = cell(size(expname));
correcttrialind = cell(size(expname));
falsetrialind = cell(size(expname));
misstrialind = cell(size(expname));
gotrial = cell(size(expname));
nogotrial = cell(size(expname));
alltrialgn = cell(size(expname));
gotrialind = cell(size(expname));
nogotrialind = cell(size(expname));
for x=1:1:length(expname)
    [kwe_infos, messages, ttl] = read_kwe_file(path2kwe{x});
    SSmessages{x} = messages(messages.nodeID == Netevent_nodeID,:);
    if cut_rec
        SSmessages{x} = SSmessages{x}(((SSmessages{x}.time_samples - start_time{x}) <= tEnd_cut) ,:);
    end
    if length(find(cellfun(@(x) numel(x) , strfind(SSmessages{x}(:,:).text,'visual')))) > 0 %visual task
        correcttrial{x}=find(strncmp(SSmessages{x}(:,:).text,'New command:SS FORWARD visual trial: correct',length('New command:SS FORWARD visual trial: correct')));
        falsetrial{x}=find(strncmp(SSmessages{x}(:,:).text,'New command:SS FORWARD visual trial: false alarm',length('New command:SS FORWARD visual trial: false alarm')));
        misstrial{x}=find(strncmp(SSmessages{x}(:,:).text,'New command:SS FORWARD visual trial: miss',length('New command:SS FORWARD visual trial: miss')));
        alltrialcf{x}=find(strncmp(SSmessages{x}(:,:).text,'New command:SS FORWARD visual trial: ',length('New command:SS FORWARD visual trial: ')));
        correcttrialind{x}=arrayfun(@(y) find(alltrialcf{x}==y), correcttrial{x});
        falsetrialind{x}=arrayfun(@(y) find(alltrialcf{x}==y), falsetrial{x});
        misstrialind{x}=arrayfun(@(y) find(alltrialcf{x}==y), misstrial{x});
        
        gotrial{x}=find(strncmp(SSmessages{x}(:,:).text,'SST PRESENT 13',length('SST PRESENT 13')));
        nogotrial{x}=find(strncmp(SSmessages{x}(:,:).text,'SST PRESENT 15',length('SST PRESENT 15')));
        alltrialgn{x}=find(strncmp(SSmessages{x}(:,:).text,'SST PRESENT 1',length('SST PRESENT 1')));
        gotrialind{x}=arrayfun(@(y) find(alltrialgn{x}==y), gotrial{x});
        nogotrialind{x}=arrayfun(@(y) find(alltrialgn{x}==y), nogotrial{x});
    
        
    elseif length(find(cellfun(@(x) numel(x) , strfind(SSmessages{x}(:,:).text,'auditory')))) > 0 %audiory task
        
         gotrial{x}=find(strncmp(SSmessages{x}(:,:).text,'SST PRESENT 13',length('SST PRESENT 13')));
        nogotrial{x}=find(strncmp(SSmessages{x}(:,:).text,'SST PRESENT 15',length('SST PRESENT 15')));
        alltrialgn{x}=find(strncmp(SSmessages{x}(:,:).text,'SST PRESENT 1',length('SST PRESENT 1')));
        gotrialind{x}=arrayfun(@(y) find(alltrialgn{x}==y), gotrial{x});
        nogotrialind{x}=arrayfun(@(y) find(alltrialgn{x}==y), nogotrial{x});
    else
        error('weird task?')
    end
end
