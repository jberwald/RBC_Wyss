% The subroutine used to condense the movie
% written by Jia-Rong Yeh  January 2009.
function [new_bund,new_img]=cond_img(img_bund,img)
img_siz=size(img_bund);
cnt=0;
new_bund=zeros(floor(img_siz(1)/3),floor(img_siz(2)/3),'int16');
for i=1:floor(img_siz(1)/3)
    for j=1:floor(img_siz(2)/3)
        s=zeros(9,5000);
        if img_bund(i*3-2,j*3-2) > 0
           s(1,:)=img(:,img_bund(i*3-2,j*3-2))';
        end
        if img_bund(i*3-2,j*3-1) > 0
           s(2,:)=img(:,img_bund(i*3-2,j*3-1))';
        end
        if img_bund(i*3-2,j*3) > 0
           s(3,:)=img(:,img_bund(i*3-2,j*3))';
        end
        if img_bund(i*3-1,j*3-2) > 0
           s(4,:)=img(:,img_bund(i*3-1,j*3-2))';
        end
        if img_bund(i*3-1,j*3-1) > 0
            s(5,:)=img(:,img_bund(i*3-1,j*3-1))';
        end
        if img_bund(i*3-1,j*3) > 0
            s(6,:)=img(:,img_bund(i*3-1,j*3))';
        end
        if img_bund(i*3,j*3-2) > 0
            s(7,:)=img(:,img_bund(i*3,j*3-2))';
        end
        if img_bund(i*3,j*3-1) > 0
            s(8,:)=img(:,img_bund(i*3,j*3-1))';
        end
        if img_bund(i*3,j*3) > 0
            s(9,:)=img(:,img_bund(i*3,j*3))';
        end
        B=s > 0;
        B_cnt=sum(B);
        if sum(B_cnt) > 30000
            cnt=cnt+1;
            new_bund(i,j)=cnt;
            new_img(:,cnt)=(sum(s)./B_cnt)';
            if mod(cnt,10) == 0
               fprintf('-');
            end
            if mod(cnt,500) ==0
                fprintf('\n');
            end
        end
        clear s B B_cnt;
    end
end