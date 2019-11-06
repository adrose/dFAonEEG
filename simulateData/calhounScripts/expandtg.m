function [tg tag] = expandtg(tg);
%% this function expands a structure tg.lab{} and tg.val{} creating
%% varibles (e.g. tg.lab{1}) with vectors (e.g. tg.val{1}). This makes it
%% easier to carry around variables and create unique folder/file
%% identifiers for the data generated with the particular variables

nvar = length(tg.val);
for G = 1:nvar
numvals(G)=length(tg.val{G});
end;

numvals = [1 numvals];
nval = prod(numvals);

cnt=0;
for G = 1:nvar
if numvals(G+1)>1;cnt=cnt+1;end;
siz1 = prod(numvals(1:G));
tg.val{G} = reshape(repmat(tg.val{G},siz1,nval/siz1/numvals(G+1)),1,nval);
end;

%% create a tag for files/folders
for G = 1:length(tg.val{1})
temptag = [];
for k = 1:length(tg.val)
    
    if iscell(tg.val{k}(G))==0
    ttag = sprintf('tg%d_',tg.val{k}(G));  
    elseif iscell(tg.val{k}(G))==1
    ttag = sprintf('tgD_'); %% D as in default parameter 
    else;end;
    
temptag = [temptag ttag];
end;
tag{G}=temptag;
end;