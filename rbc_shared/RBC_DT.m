% This program is written by Jia-Rong Yeh on someday of January, 2009
% It transfers the data of 10 new red blood cells (RBC) saved in subdirectories
% named from new_1 to new_10. Each subdirectory includes 5000 frame from
% 0000 to 4999, which was saved as new_?0000.txt to new_?4999.
% In a frame, the values of pixels larger than 0 mean which are inside the
% boundary of RBC.

clear;
for i=1:1 % 10   % set the number of new RBC
    subdir=sprintf('new_%d',i); % assemble the name of subdirectory
    cd(subdir); % change current directory to the subdirectory of cell
    % ==========  assemble the file name of frame 
    for j=0:4999  
        if j < 10
            filnam=sprintf('new_%d000%d.txt',i,j);
        else
            if j < 100
                filnam=sprintf('new_%d00%d.txt',i,j);
            else
                if j < 1000
                    filnam=sprintf('new_%d0%d.txt',i,j);
                else
                    filnam=sprintf('new_%d%d.txt',i,j);
                end
            end
        end
        % ====================== open the txt file of frame
        fid=fopen(filnam,'r');
        if fid > 0   % check the existence of frame 
            fprintf('Now is opening %s => ',filnam);
            img_siz=fscanf(fid,'%d %d');   % read image size
            c=fread(fid,2);
            tmp_v=fscanf(fid,'%d');
            if j==0   % ====================== find the boundary of cell when frame 1
                cnt1=0; cnt2=0;
                img_bund=zeros(img_siz(2),img_siz(1)); % allocate the memory of boundary matrix
                for n=1:img_siz(2)
                    for m=1:img_siz(1)
                        cnt1=cnt1+1;
                        if tmp_v(cnt1) >0
                            cnt2=cnt2+1;  % count the number of pixels inside the cell
                            img_bund(m,n)=cnt2; % give a ID number to the pixel inside the cell
                        end
                    end
                end
                
                img=zeros(5000,cnt2,'int16'); % allocate the memory for 5000 frames
                % ======================== map frame 1 to img(1,:) using
                % ID of pixels
                cnt1=0;
                for m=1:img_siz(2)
                    for n=1:img_siz(1)
                        cnt1=cnt1+1;
                        if img_bund(m,n) > 0
                            img(j+1,img_bund(m,n))=tmp_v(cnt1);
                            C(m,n)=tmp_v(cnt1);
                        end
                    end
                end
            else
                % map frame 2 to 50000 to img(2:5000,:)
                cnt1=0;
                for n=1:img_siz(2)
                    for m=1:img_siz(1)
                        cnt1=cnt1+1;
                        if img_bund(m,n) > 0
                            img(j+1,img_bund(m,n))=tmp_v(cnt1);
                        end
                    end
                end
            end
            %fprintf('The image size is %d * %d . string is %s \n',img_siz,c);
            fclose(fid);
        end
        fprintf('The frame %d is completed.\n',j);
    end
    % ============== show the frame1 of cell 
    figure(1);
    image(C/60);colorbar;
    cd ..;  % change the current directory to the original directory
    filnam=sprintf('new_%d.mat',i);  % assemble the file name of *.mat
    save(filnam,'img','img_bund');   % save the data
end