
%indicate whether xml file is okay and if not whether BOT was used

prompt = {'xml file okay (0) or corrupt (1) ?'};
dlg_title = '0 for xml okay, 1 for corrupt';
num_lines = 1;
def = {'0'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
answer = str2mat(answer);
xml_file_corrupt = str2num(answer);


if xml_file_corrupt
prompt = {'Used BOT (1) or not (0) ?'};
dlg_title = '1 if BOT, 0 if not';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
answer = str2mat(answer);
BOT = str2num(answer);
end

[xml_FileName,PathName] = uigetfile('*.xml','Select the xml-file');

addpath(PathName);
Loc  = strcat(PathName,xml_FileName);
Folder = findstr(Loc,'\');
Folder = Loc(1:Folder(end));



if xml_file_corrupt
xml_FileName = dir(strcat(Folder,'\*C*.xml'));
xml_FileName_aux = dir(strcat(Folder,'\x*.xml'));
else
    xml_FileName = dir(strcat(Folder,'\*.xml'));
end


%read xml-file(s)

if xml_file_corrupt
    xml_file_aux = xml_parseany(fileread(xml_FileName_aux.name));
    
    if BOT
        xml_file = xml_parseany(fileread(xml_FileName.name));
        initial.max_frame = str2num(xml_file.TSeries{1,1}.PVTSeriesElementSequenceBOT{1,1}.ATTRIBUTE.repetitions);
        initial.numberofframes = initial.max_frame;
        initial.msperline = 1000*str2num(xml_file.TSeries{1,1}.PVTSeriesElementSequenceBOT{1,1}.ATTRIBUTE.repetitionPeriod);
        
    else
        xml_file = xml_parseany(fileread(xml_FileName.name));
        initial.max_frame = str2num(xml_file.TSeries{1,1}.PVTSeriesElementSequenceTimed{1,1}.ATTRIBUTE.repetitions);
        initial.numberofframes = initial.max_frame;
        initial.msperline = 1000*str2num(xml_file.TSeries{1,1}.PVTSeriesElementSequenceTimed{1,1}.ATTRIBUTE.repetitionPeriod);
    end
    
else
    xml_file = xml_parseany(fileread(xml_FileName.name));
    initial.max_frame = max(size(xml_file.Sequence{1,1}.Frame));
    initial.numberofframes = initial.max_frame;
    initial.msperline = -1000*(str2num(xml_file.Sequence{1,1}.Frame{1,1}.ATTRIBUTE.relativeTime) - str2num(xml_file.Sequence{1,1}.Frame{1,initial.max_frame}.ATTRIBUTE.relativeTime))/...
        (initial.max_frame-1);
end



%check for slow frames


if xml_file_corrupt ==0
    for i = 1:size(initial.tiffmovie,4)
        initial.frame_times(i) = str2num(xml_file.Sequence{1,1}.Frame{1,i}.ATTRIBUTE.relativeTime);
    end
end

frame_times_diff = diff(diff(initial.frame_times));

results.slow_frames = 0;
results.slow_frame_vector = 0;
results.number_slow_frames = 0;
x = 1;
if xml_file_corrupt ==0
    for i = 1:size(initial.tiffmovie,4)-2
        if frame_times_diff(i)>0.00001
            results.slow_frame_vector(i) = 1;
            results.slow_frames(x) = i;
            results.number_slow_frames = results.number_slow_frames + 1;
            x = x+1;
        end
    end
end



%record initial settings from xml-file(s) (not finished for the case where xml file is corrupt) 


if xml_file_corrupt
initial.date_time = xml_file_aux.ATTRIBUTE.date;
initial.microns_per_pixel = str2num(xml_file.PVStateShard{1,1}.Key{1,64}.ATTRIBUTE.value);
initial.magnification = str2num(xml_file.PVStateShard{1,1}.Key{1,187}.ATTRIBUTE.value);
initial.zoom = str2num(xml_file.PVStateShard{1,1}.Key{1,11}.ATTRIBUTE.value);
   
else
initial.date_time = xml_file.ATTRIBUTE.date;
initial.microns_per_pixel = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,22}.ATTRIBUTE.value);
initial.magnification = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,6}.ATTRIBUTE.value);
initial.zoom = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,21}.ATTRIBUTE.value);
initial.numerical_aperture = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,5}.ATTRIBUTE.value);
initial.dwell_time_us = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,14}.ATTRIBUTE.value);
initial.pockels_power = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,28}.ATTRIBUTE.value);
initial.x_position = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,16}.ATTRIBUTE.value);
initial.y_position = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,17}.ATTRIBUTE.value);
initial.z_position = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,18}.ATTRIBUTE.value);
initial.PMT_voltage_1 = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,24}.ATTRIBUTE.value);
initial.PMT_voltage_2 = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,25}.ATTRIBUTE.value);
initial.PMT_preamp_gain_1 = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,32}.ATTRIBUTE.value);
initial.PMT_preamp_gain_2 = str2num(xml_file.Sequence{1,1}.Frame{1,1}.PVStateShard{1,1}.Key{1,33}.ATTRIBUTE.value);
end