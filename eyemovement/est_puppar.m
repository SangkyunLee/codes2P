function [fitInfo, fail] = est_puppar(data, CM,sPIX,pupil_radius,htwin,frames,mode)

if nargin<5,
    htwin=0;
end

Nframe = size(data,3);
if nargin<6,
    frames = 1:Nframe;
end
if nargin<7,
    mode = 1;
end
s = struct('par',[],'data',[]);
fitInfo = repmat(s,Nframe,1);
error_list = zeros(Nframe,1);

Nlen = length(frames);
fitInfo0 = repmat(s,Nlen,1);
error_list0 = zeros(Nlen,1);
parfor iimg0 = 1:Nlen   
    data1=data;
    iimg = frames(iimg0);
    
        
%     if iimg==445
%         keyboard;
%     end
    
    CMi = CM(:,iimg);
    if (iimg-htwin)<1 || (iimg+htwin)>Nframe        
        A = data1(:,:,iimg);
    else
        A = mean(data1(:,:,iimg-htwin:iimg+htwin),3);
    end
    A1 = 255- A;


    center=[round(CMi(2)) round(CMi(1))];

    try
       

        if mode==1
            pixel_list = search_pupbound(A1,center,pupil_radius,'threshold');       
            
            msk = zeros(size(A1));            
            pixel_list = union(pixel_list,sPIX{iimg});               
            msk(pixel_list)=1;    
            [X0, Y0]=cmass(msk);                    
            pixel_list = search_pupbound(A1,[Y0, X0],pupil_radius,'threshold');
            
        elseif mode==2 % do not rely the sPIX
            pixel_list = search_pupbound(A1,center,pupil_radius,'threshold');                   
        elseif mode==3 % do not rely search_pupbound
            pixel_list = sPIX{iimg};   
        end

        msk = zeros(size(A1));
        msk(pixel_list)=1;

        bw=bwperim(msk,4);
        boundary=find(bw(:));

        [y, x]= ind2sub(size(A1), boundary);
        if mode==3,
            k = convhull(x,y);
            x = x(k);
            y = y(k);
            %figure; plot(x,y,'r-')
            boundary = sub2ind(size(A1),y,x);
        end

        [z, r] = fitcircle([x y]');
        fitInfo0(iimg0).par=[z;r];
        fitInfo0(iimg0).data = boundary;
    catch
        error_list0(iimg0)=1;
    end
    
end
fitInfo(frames) = fitInfo0;
error_list(frames) = error_list0;


fail = error_list;