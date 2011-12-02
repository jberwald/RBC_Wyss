% this program is used to decompose the first 7 mode function of a RBC via the condensed data. 
% written by Jia-Rong Yeh on February 2009.

clear;

% need to update file names, number of frames (5000); sum (B)

filnam=sprintf('C_4_cond.mat');
load(filnam);
img_siz=size(new_bund);
mem_scale=size(new_img);
point_num=mem_scale(2);
MD=zeros(5000,point_num,7,'single');
for i=1:img_siz(1)
    for j=1:img_siz(2)
        if new_bund(i,j) > 0
            B=new_img(1:5000,new_bund(i,j)) > 0;
            if sum(B) > 4000 %smaller than the number of frames
                fprintf('.'); 
                [c,r]=emd(1:5000,new_img(1:5000,new_bund(i,j))',5);
                for k=1:5
                    MD(:,new_bund(i,j),k)=single(c(k,:)');
                end
                clear c r;
                if mod(new_bund(i,j),50) == 0
                    fprintf(' %d\n',new_bund(i,j));
                end
            end
        end
    end
end
B=isinf(MD);
MD(B)=0;
filnam=sprintf('C_4_MF.mat');
save(filnam,'new_bund','MD');
clear new_bund MD;