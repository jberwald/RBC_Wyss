clear;
new_id=input('Which sample of new RBC do you want? (from 1 to 10 excluding 3 and 8):');
old_id=input('Which sample of old RBC do you want? (from 1 to 10): ');

% ====================load new RBC =========================
filnam=sprintf('A_4_MF.mat',new_id);
fprintf(' Complement level: High-----------------------\n');
load(filnam);
pix1_cnt=max(max(new_bund));
%MD=double(MD);
for i=1:4
    Max(i)=max(max(MD(:,:,i)));
    Min(i)=min(min(MD(:,:,i)));
    fprintf('For mode %d, max is %f and min is %f.\n',i,Max(i),Min(i));
    MD(:,:,i)=(MD(:,:,i)-Min(i))/(Max(i)-Min(i));
    %mean_v=mean(mean(MD(:,:,i)));
    %std_v=std(std(MD(:,:,i)));
    MD(:,:,i)=64*MD(:,:,i);
end
MF1=MD(:,:,1:4);
bound1=new_bund;
clear MD new_bund;
% ====================load old RBC =========================
filnam=sprintf('C_4_MF.mat',old_id);
fprintf(' Complement level: Low-----------------------\n');
load(filnam);
pix2_cnt=max(max(new_bund));
%MD=double(MD);
for i=1:4
    Max(i)=max(max(MD(:,:,i)));
    Min(i)=min(min(MD(:,:,i)));
    fprintf('For mode %d, max is %f and min is %f.\n',i,Max(i),Min(i));
    MD(:,:,i)=(MD(:,:,i)-Min(i))/(Max(i)-Min(i));
    %mean_v=mean(mean(MD(:,:,i)));
    %std_v=std(std(MD(:,:,i)));
    MD(:,:,i)=64*MD(:,:,i);
end
MF2=MD(:,:,1:4);
bound2=new_bund;
clear MD new_bund;
fprintf('The mode functions have been loaded.\n');

% Preprosessing --------------------------------
n1=max(max(bound1));
px1=zeros(1,n1,'int16');
py1=zeros(1,n1,'int16');
img_siz1=size(bound1);
for i=1:img_siz1(1)
    for j=1:img_siz1(2)
        if bound1(i,j) > 0
            px1(bound1(i,j))=i;
            py1(bound1(i,j))=j;
        end
    end
end
C1=zeros(img_siz1(1),img_siz1(2),4);

n2=max(max(bound2));
px2=zeros(1,n2,'int16');
py2=zeros(1,n2,'int16');
img_siz2=size(bound2);
for i=1:img_siz2(1)
    for j=1:img_siz2(2)
        if bound2(i,j) > 0
            px2(bound2(i,j))=i;
            py2(bound2(i,j))=j;
        end
    end
end
C2=zeros(img_siz2(1),img_siz2(2),4);


% std_def=0.02*ones(1,7);
for i=1:600
    fig=figure(1);
    for j=1:4
        mean_C=mean(MF1(i,:,j));
        std_C=std(MF1(i,:,j));
        for k=1:n1
            C1(px1(k),py1(k),j)=16*((MF1(i,k,j)-mean_C)/std_C)+32;
        end
        txt=sprintf('Mode %d',j);
        subplot(2,4,j);
        image(C1(:,:,j));
        xlabel(num2str(std_C));
        title(txt);
        if j==1
            ylabel('Pathologic RBC');
        end
    end
    for j=1:4
        mean_C=mean(MF2(i,:,j));
        std_C=std(MF2(i,:,j));
        for k=1:n2
            C2(px2(k),py2(k),j)=16*((MF2(i,k,j)-mean_C)/std_C)+32;
        end
        txt=sprintf('Mode %d',j);
        subplot(2,4,j+4);
        image(C2(:,:,j));
        xlabel(num2str(std_C));
        title(txt);
        if j==1
            ylabel('Healthy RBC');
        end
    end
     F(i)=getframe(fig);
end
movie2avi(F,'RBC.avi','compression','Cinepak');