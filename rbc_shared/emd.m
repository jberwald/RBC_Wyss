function [c,r]=emd(s,r,nlayer)
n=length(r);
% ======================= decomposing process==================
for i=1:nlayer
    SD=1;
    h=r;
    rest=1;
    for  k=1:7 % new criterion 
        % ========================== find down_envelope ==========
        cnt=1;
        clear extrema x;
        extrema(cnt)=h(1);
        x(cnt)=1;
        for j=2:n-1
            if h(j) < h(j-1) & h(j) < h(j+1)
                cnt=cnt+1;
                extrema(cnt)=h(j);
                x(cnt)=j;
            end
        end
        cnt=cnt+1;
        extrema(cnt)=h(j);
        x(cnt)=n;
        pp=spline(x,extrema);
        dv=ppval(pp,s);
        % ========================== find up_envelope ==========
        cnt=1;
        clear extrema x;
        extrema(cnt)=h(1);
        x(cnt)=1;
        for j=2:n-1
            if h(j) > h(j-1) & h(j) > h(j+1)
                cnt=cnt+1;
                extrema(cnt)=h(j);
                x(cnt)=j;
            end
        end
        cnt=cnt+1;
        extrema(cnt)=h(j);
        x(cnt)=n;
        pp=spline(x,extrema);
        s=1:n;
        uv=ppval(pp,s);
        m=(uv+dv)/2;
        
        %fprintf('h %d\n', size(h));
        
        h=h-m;
        rest=0;
    end
    c(i,:)=h;
    r=r-h;
end