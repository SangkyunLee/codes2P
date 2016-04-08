function [pixel_list,boundary_list] = search_pupbound(map,center,search_radius,method)
% function [pixel_list,boundary_list] = search_pupbound(map,center,search_radius,method)
% This function is a little modified version of disk_roi.m of Tsai-Wen Chen
%    

%

ntheta=90;
nsample = search_radius;
t=(1:ntheta)/ntheta*2*pi;


ss=size(map);
map=imresize(map,ss);
map=double(map);


x=round(center(2));
y=round(center(1));

%---- est init search bound
lprofile0=zeros(nsample,ntheta);
for i=1:ntheta
    f=improfile(map,[x,x+search_radius*cos(t(i))],[y,y+search_radius*sin(t(i))],nsample,'bilinear')';    
    lprofile0(:,i)=f;
end
mL= mean(lprofile0,2);
%[~,full_radius] = min(diff(mL));
full_radius= find_first_cross(mL,0.5);   

full_radius = round(full_radius*1.5);



lprofile=zeros(nsample,ntheta);
lprofile_mod=zeros(nsample,ntheta);
r_threshold_cross=zeros(1,ntheta);
r_stiff=zeros(1,ntheta);
for i=1:ntheta
    f=improfile(map,[x,x+full_radius*cos(t(i))],[y,y+full_radius*sin(t(i))],nsample,'bilinear')';    
    lprofile(:,i)=f;

    
    switch method
        case 'threshold'
            r_threshold_cross(i)=find_first_cross(f,0.8);       
            lprofile_mod(:,i)=-1*abs((1:nsample)-r_threshold_cross(i))+nsample;        
        case 'slope'
            [~,inx] = min(diff(f));
            r_stiff(i)= inx;
            lprofile_mod(:,i)=-1*abs((1:nsample)-r_stiff(i))+nsample;
    end
end

path=find_path(lprofile_mod);
rin=round(path);
%% determine pixels

roimap=zeros(ss);
for xx=-full_radius:full_radius
    for yy=-full_radius:full_radius
        
        theta=angle(xx+yy*sqrt(-1));
        if theta<=0;
            theta=theta+2*pi;
        end
        rr=sqrt(xx^2+yy^2)/full_radius*nsample;
        tind=ceil(theta/(2*pi)*ntheta);
        rin_int=rin(tind);
       

        indx=xx+x;
        indy=yy+y;
        if (indx>0) && (indy>0) && indx<=(ss(2))  && indy<=(ss(1))
            if rr<rin_int
                roimap(indy,indx)=1;
            end
        end
    end
end



pixel_list=find(roimap(:)==1);
       
bw=bwperim(roimap,4);
boundary_list=find(bw(:));




end




function r_threshold_cross=find_first_cross(f,th, minf)
    if nargin<3
        globalmin=min(f);
    else
        globalmin =minf;
    end
    
    threshold=(max(f(1:5))-globalmin)*th+globalmin;
    r_threshold_cross=find(f<threshold,1,'first')-1;

end
    


function path=find_path(lprofile)
    pointer=zeros(size(lprofile));
    value=zeros(size(lprofile));
    value(:,1)=lprofile(:,1);
    for i=2:size(lprofile,2)
        for j=2:(size(lprofile,1)-1)
            [M,ind]=max(value((j-1):(j+1),i-1));
            value(j,i)=M+lprofile(j,i);
            pointer(j,i)=j+ind-2;
        end
    end

    %%second traverse to minimize boundary effect
    pointer=zeros(size(lprofile));
    value(:,1)=value(:,end);
    for i=2:size(lprofile,2)
        for j=2:(size(lprofile,1)-1)
            [M,ind]=max(value((j-1):(j+1),i-1));
            value(j,i)=M+lprofile(j,i);
            pointer(j,i)=j+ind-2;
        end
    end
    
    path=zeros(1,size(lprofile,2));
    [M,ind]=max(value(:,end));
    path(end)=ind;
    
    for j=(size(lprofile,2)):-1:2
    path(j-1)=pointer(path(j),j);
    end
end