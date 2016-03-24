function s=nex2mat(Folder,nex_file)
%retrieve trial info from *.nex file in matlab
%
%nex file
[nex_folder,nex_filename,nex_ext]=fileparts(nex_file);
nex_filename = [nex_filename nex_ext];
s=struct;
s.nexFile = nex_filename;
s.nexFolder = nex_folder;
s.nexData = readNexFile(fullfile(Folder,s.nexFile));


