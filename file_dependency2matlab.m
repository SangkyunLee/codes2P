
fnpath = '/home/slee/data/codes2P/cellidentification/roigui/roigui(release)/private';
%     'Image Processing Toolbox'    'MATLAB'
fnpath ='/home/slee/data/codes2P/cellidentification'
% 'Image Processing Toolbox'    'MATLAB'    'Signal Processing Tool…' 

fnpath ='/home/slee/data/codes2P/eyemovement'
%     'Computer Vision System Toolbox'
%     'Control System Toolbox'
%     'Image Processing Toolbox'
%     'MATLAB'
%     'Statistics and Machine Learning Toolbox'
fnpath ='/home/slee/data/codes2P/analysis'
%     'MATLAB'
%     'Statistics and Machine Learning Toolbox'
fnpath ='/home/slee/data/codes2P/deconvolution'
% 'MATLAB'
fnpath ='/home/slee/data/codes2P/extras'
fnpath ='/home/slee/data/codes2P/preprocessing'
%     'Image Processing Toolbox'
%     'MATLAB'
%     'Signal Processing Toolbox'
fnpath='/home/slee/data/codes2P/smlr'
fnpath='/home/slee/data/codes2P/efficient_subpixel_registration'
%     'MATLAB'
% when I use gpu, following toolbox requires 
% 'Image Processing Toolbox'        'Parallel Computing Too…' 
fnpath='/home/slee/data/codes2P/plots'

fnpath ='/home/slee/data/codes2P/'

flist = dir([fnpath filesep '*.m']);
nf =length(flist);
files_in = cell(nf,1);
for i = 1 : nf
files_in{i}=fullfile(fnpath,flist(i).name);
end

[names, folders] = dependencies.toolboxDependencyAnalysis(files_in(1)); 

namessub= cell(length(files_in),1);
inx = [];
for i = 1: length(files_in)
    
    namessub{i} = dependencies.toolboxDependencyAnalysis(files_in(i));
    if length(namessub{i})==1,
        namessub{i} =[];
    else
        inx=[inx i];
    end
end
7, 9
%% List of toolboxes needed
%'Image Processing Toolbox'
%'Statistics and Machine Learning T??
%'Signal Processing Toolbox'
%'Computer Vision System ??
%'parallel computing toolbox'
% gpu_LDA.m need control toolbox because it uses parallel.gpu.


% for eyetracking
% I need Data Acquisition Toolbox and Imagge Acquisition Toolbox
