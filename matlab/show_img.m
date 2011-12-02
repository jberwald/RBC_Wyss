% This program is used to play the movies of two red blood cells (RBC).
% Data of movies are saved as new_?.mat or old_?.mat.
% Author: Jia-Rong Yeh  January 17, 2009

clear;
id1=0;  % asign the sample ID of new and old RBC
while id1 < 1 | id1 > 10 
     id1=input('Which sample of New RBC you want to watch? (from 1 to 10) ');
end 
id2=0;
while id2 < 1 | id2 > 10
    id2=input('Which sample of Old RBC you want to watch? (from 1 to 10) ');
end

% ============== associate the file name of data ====================
filnam=sprintf('new_cell.mat');  % the file name contains 'norm' means it is a normalized data
% In a normalized data, the value of data locates on range from 0 to 1.
% You can also display the original (not normalized) data by removing the
% characters 'norm'
load(filnam);
img_siz1=size(img_bund);
img_bund1=img_bund;
img1=img;
filnam=sprintf('old_%dnorm.mat',id2);
load(filnam);
img_siz2=size(img_bund);

% ======================= play movie ===========================
for i=1:5000    % i is the frame id
    C1=zeros(img_siz1(1),img_siz1(2),'int16');
    C=zeros(img_siz2(1),img_siz2(2),'int16');
    % ==============================  revise the row vector to a frame
    for j=1:img_siz1(1)
        for k=1:img_siz1(2)
            if img_bund1(j,k) > 0
                C1(j,k)=img1(i,img_bund1(j,k))*128;
            end
        end
    end
    
    for j=1:img_siz2(1)
        for k=1:img_siz2(2)
            if img_bund(j,k) > 0
                C(j,k)=img(i,img_bund(j,k))*128;
            end
        end
    end
    % ========== show the frames ==============
    figure(1);
    fram_txt=sprintf('Frame %d',i);
    subplot(121);
    title(fram_txt);
    image(C1);
    txt=sprintf('New RBC %d',id1);
    xlabel(txt);
    subplot(122);
    title(fram_txt);
    image(C)
    txt=sprintf('Old RBC %d',id2);
    xlabel(txt);
end