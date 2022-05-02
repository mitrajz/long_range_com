
cd(['/mnt/data/Mitra/raw/',animalname,'/oemessages'])

SSmessages = readtable(oemessagesfilename,'Delimiter',',','HeaderLines', 0);
SSmessages = table2cell(SSmessages); % all SS messages
if ~rotate
    % index of the start of 10-28 part
    startindex = find(strncmp(SSmessages,'Go stim id is',length('Go stim id is')))
    
    
    if strcmp(SSmessages(startindex),'Go stim id is 28')
        gostimid = '28';
        nogostimid = '10';
    elseif strcmp(SSmessages(startindex),'Go stim id is 10')
        gostimid = '10';
        nogostimid = '28';
    end
    % make sure no change of go/nogo stim id after this point
    %%%
    gotrial=find(strncmp(SSmessages,['SS PRESENT ',gostimid,'.'],length(['SS PRESENT ',gostimid,'.'])));
    nogotrial=find(strncmp(SSmessages,['SS PRESENT ',nogostimid,'.'],length(['SS PRESENT ',nogostimid,'.'])));
    alltrialgn=find(strncmp(SSmessages,'SS PRESENT',length('SS PRESENT')));
    gotrialind=arrayfun(@(y) find(alltrialgn==y), gotrial);
    nogotrialind=arrayfun(@(y) find(alltrialgn==y), nogotrial);
else
    % index of the start of 1-19 part
    startindexp = find(strncmp(SSmessages,'Go stim id changed to 1',length('Go stim id changed to 1')));
    startindexn = find(strncmp(SSmessages,'Go stim id changed to 19',length('Go stim id changed to 19')));
    startindex = [startindexp startindexn]
    
    if strcmp(SSmessages(startindex),'Go stim id changed to 1')
        gostimid = '1';
        nogostimid = '19';
    elseif strcmp(SSmessages(startindex),'Go stim id changed to 19')
        gostimid = '19';
        nogostimid = '1';
    end
    % make sure no change of go/nogo stim id after this point
    %%%
    gotrial=find(strncmp(SSmessages,['SS PRESENT ',gostimid,'.'],length(['SS PRESENT ',gostimid,'.'])));
    nogotrial=find(strncmp(SSmessages,['SS PRESENT ',nogostimid,'.'],length(['SS PRESENT ',nogostimid,'.'])));
    alltrialgn=find(strncmp(SSmessages,'SS PRESENT',length('SS PRESENT')));
    gotrialind=arrayfun(@(y) find(alltrialgn==y), gotrial);
    nogotrialind=arrayfun(@(y) find(alltrialgn==y), nogotrial);
    
    
end
