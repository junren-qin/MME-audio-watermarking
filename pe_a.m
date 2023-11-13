function m=pe_a(B,R,alpha,sita_n,type)
    N=size(B,1);
    r1=eqvalent_radius((1-alpha)*(1/type)*R*B);
    r2=eqvalent_radius(B);
    syms n x;
    c1=(r1^2-r2^2+sita_n*n)/(2*r1*sqrt(sita_n*n));
    c2=(r2^2-r1^2+sita_n*n)/(2*r2*sqrt(sita_n*n));
    h1=r1*(1-c1);
    h2=r2*(1-c2);
    funb=(cos(x))^N;
    l1=asin(1-h1/r1);
    l2=asin(1-h2/r2);
    beta1=int(funb,x,l1,pi/2);
    beta2=int(funb,x,l2,pi/2);
    V1=sphere_volume(r1,N)-sphere_parameter(r1,N)*beta1+sphere_parameter(r2,N)*beta2;
    V2=sphere_parameter(r1,N)*beta1+sphere_parameter(r2,N)*beta2;
    fun=chi2_pdf(N,n);
    fun3=matlabFunction(V1*fun);
    fun4=matlabFunction(V2*fun);
    AV=integral(matlabFunction(fun),0,(r2-r1)^2/sita_n)*sphere_volume(r1,N)...
        +integral(fun3,(r2-r1)^2/sita_n,(r2^2-r1^2)/sita_n)...
        +integral(fun4,(r2^2-r1^2)/sita_n,(r1+r2)^2/sita_n); 
%     AV=quad(fun3,(r2-r1)^2/sita_n,(r1+r2)^2/sita_n);
    m=AV/det((1-alpha)*(1/type)*R*B);
end

function m=eqvalent_radius(B)
    N=size(B,1);
    m=pi^(-0.5)*gamma(1+N*.5)^(1/N)*det(B)^(1/N);
end

function m=inc_betafunction(h,N)
    syms x;
    fun=x^((N-1)/2)*(1-x)^(-.5);
    m=int(fun,x,0,h);
end

function m=sphere_parameter(r,N)
    m=r^N*pi^((N-1)*.5)/(gamma((N+1)*.5));
end

function m=chi2_pdf(N,x)
    m=((.5)^(N/2)/gamma(N/2))*x^(N/2-1)*exp(1)^(-x/2);
end

function m=sphere_volume(r,N)
    m=r^N*pi^((N)*.5)/(gamma((N+2)*.5));
end

