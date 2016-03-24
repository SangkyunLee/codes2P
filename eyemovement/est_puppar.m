function [fitInfo, fail] = est_puppar(data, CM,sPIX)

Nframe = size(data,3);
s = struct('par',[],'data',[]);
fitInfo = repmat(s,Nframe,1);
pupil_radius=100;

error_list = zeros(Nframe,1);
parfor iimg = 1:Nframe   
    CMi = CM(:,iimg);
    A = data(:,:,iimg);
    A1 = 255- A;


    center=[round(CMi(2)) round(CMi(1))];
    
    try
        pixel_list = disk_roi(A1,center,pupil_radius,[1 1],0);
        pixel_list = union(pixel_list,sPIX{iimg});   
    
        msk=zeros(size(A));
        msk(pixel_list)=1;    
        [X0, Y0]=cmass(msk);
        
    
        pixel_list = disk_roi(A1,[Y0, X0],pupil_radius,[1 1],0);
        pixel_list =union(pixel_list, sPIX{iimg});    

        msk = zeros(size(A1));
        msk(pixel_list)=1;

        bw=bwperim(msk,4);
        boundary=find(bw(:));

        [y, x]= ind2sub(size(A1), boundary);

        [z, r] = fitcircle([x y]');
        fitInfo(iimg).par=[z;r];
        fitInfo(iimg).data = boundary;
    catch
        error_list(iimg)=1;
    end
end
fail = error_list;