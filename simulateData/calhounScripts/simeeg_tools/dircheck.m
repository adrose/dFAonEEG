function dircheck(fname);
%% DIRCHECK checks whether a directory exists
%% If it does not exist it creates it.
%% It then cd's to that directory.
%% Input: string "fname" specifying directory name

if eval([sprintf('isdir(''%s'')',fname)])==0
eval([sprintf('mkdir %s',fname)]);
else;end;
eval([sprintf('cd %s',fname)]);