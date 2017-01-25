function mp = get_OSpath(mp,mp_os)
% function mp = get_OSpath(mp,mp_os)
% get a new OS dependent path.
% mp_os ={'/media/sdb_WD4T/data_2photon','Y:/data_2photon'}
% The order of OS is unix_path, win_path....


if strfind(computer,'PC') && ~isempty(strfind(mp,mp_os{1}))
    mp = strrep(mp,mp_os{1},mp_os{2});
elseif strfind(computer,'GLNX') && ~isempty(strfind(mp,mp_os{2}))
    mp = strrep(mp,mp_os{2},mp_os{1});
end