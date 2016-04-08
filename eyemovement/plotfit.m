function fiterror =plotfit(movfn, data,framelist, fitInfo)


writerObj = VideoWriter(movfn,'MPEG-4');
writerObj.FrameRate = 60;
open(writerObj);


t = linspace(0, 2*pi, 100);

fiterror = zeros(max(framelist),1);


for i = 1:length(framelist)
    iimg = framelist(i);
%     if i==358, keyboard; end
   
    P = fitInfo(iimg).par;
    A = data(:,:,iimg);            
    A1 =histeq(A);        
    boundary = fitInfo(iimg).data;
    msk = A1;
    if max(boundary)<=length(msk(:))        
        msk(boundary)=255;
    end
        
    if ~isempty(P)

        
%         if mod(iimg,5)==1.
            
%             imagesc(msk)
%             hold on;            
%             plot(P(1)  + P(3)  * cos(t), P(2)  + P(3) * sin(t), 'r')
            try 
                inx = sub2ind([size(data,1) size(data,2)],round(P(2)  + P(3)  * sin(t)), round(P(1)  + P(3)  * cos(t)));
                msk(round(inx))=255;
                position = [10 10];
                textString = sprintf('F=%d,R=%d',iimg,round(P(3)));
                msk = insertText(msk,position,textString);
            catch
                fiterror(iimg)=1;
            end
            
%             axis equal
%             title(sprintf('%d',iimg));
%             frame = getframe;
            
           
%         end
        
    end
    writeVideo(writerObj,msk);
end

close(writerObj);