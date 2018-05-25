classdef TPSession2 
    % 2018-03-02, Sangkyun Lee

    properties
        fn_data; % data_file after calling load_data
        
       
        session_path; % root_path for all session data
        
        animal_state; % 'AW','AN(fentanyl)', etc
        imagingmethod; % spiral(224), sprial(512), resonant(1024)
        
        no_scan; % # scans
        scanId; % scan identity numbers     
        scan_ncell; % no cell identified
        FOV_size; % FOV size in pixel   
        
        
        
        scans=struct([]); % data loaded by load_data function        
        rotamotion=struct([]); % motion_estimate (unit: image frame)        
        eye = struct([]); 
        vstim = struct([]);
        
    end
    methods 
        function self =TPSession2(fn_data,  info)
            if exist(fn_data,'file')~=2
               error('File not exist:%s',fn_data);
               % TO DO for later % loading 
               % I have to make a function to load data into a matlabfile
               % load_data(info)
            else
                loadscan            
            end            
        end
        
        
        
       
    end
end



    