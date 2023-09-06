function xx=chase(a,b,c,rr);
    n=length(a);
    b=[0 b];
    w(1)=rr(1)/a(1);
    v(1)=c(1)/a(1);
    for i=2:n-1
        w(i)=(rr(i)-b(i)*w(i-1))/(a(i)-b(i)*v(i-1));
        v(i)=c(i)/(a(i)-v(i-1)*b(i));
    end 
    w(n)=(rr(n)-b(n)*w(n-1))/(a(n)-b(n)*v(n-1));
    %%
    xx(n)=w(n);
    for i=n-1:-1:1
        xx(i)=w(i)-v(i)*xx(i+1);
    end