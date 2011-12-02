clear;
%for i=1:1 % 10
%subdir=sprintf('new_%d',i);
%cd(subdir);
    fid = fopen('boundary.txt');
    img_siz=fscanf(fid,'%d %d');
    c=fread(fid,2);
    tmp_v=fscanf(fid,'%d');
    cnt1=0; cnt2=0;
    img_bund=zeros(img_siz(2),img_siz(1));
    for n=1:img_siz(2)
        for m=1:img_siz(1)
            cnt1=cnt1+1;
            if tmp_v(cnt1) > 0
                cnt2=cnt2+1;
                img_bund(n,m)=cnt2;  
            end
        end
    end
    fclose(fid);
    img=zeros(4999,cnt2,'int16');
    for j=0:4999
        if j < 10
            filnam=sprintf('C_40000%d.txt',j);
        else
            if j < 100
                filnam=sprintf('C_4000%d.txt',j);
            else
                if j < 1000
                    filnam=sprintf('C_400%d.txt',j);
                else
                    filnam=sprintf('C_40%d.txt',j);
                end
            end
        end
        fid=fopen(filnam,'r');
        if fid > 0
            fprintf('Now is opening %s => ',filnam);
            img_siz=fscanf(fid,'%d %d');
            c=fread(fid,2);
            tmp_v=fscanf(fid,'%d');
            %img_bund=img_bund';
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
            %fprintf('The image size is %d * %d . string is %s \n',img_siz,c);
            fclose(fid);
            fprintf('The frame %d is completed mean value = %f .\n',j,mean(img(j+1,:)));
        end
        
    end
    figure(1);
    image(C/60);colorbar;
    % cd ..;
    filnam=sprintf('C_4.mat',i);
    save(filnam,'img','img_bund');
%end