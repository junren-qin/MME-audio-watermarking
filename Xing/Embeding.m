%%%%%%%%%%%%%%%%%%嵌水印子程序%%%%%%%%%%%%%%%%%%%%
function [ Rw,E,d ] = Embeding(Ro,wm,S,T,G,L,n)

Rw=Ro;  
Lg=length(Rw);     %采样点个数
L1=floor(Lg/S);     %L1个帧
n1=n+1;
start=1;     %嵌入起点

%%%%%%%%%%%%%%%%%%计算二阶差分统计量%%%%%%%%%%%%%%%
A=zeros(1,n1);
for i=1:n1
    A(i)=(-1)^(i+1)*nchoosek(n,i-1); %计算n阶差分方程的带符号系数
end
A1=abs(A);%计算n阶差分方程的无符号系数，嵌入过程需要用到

E=zeros(1,L);      
k=1;
for i=start:L1
    Rt=Rw((i-1)*S+1:(i-1)*S+S);
    for j=1:S/n1
        Rt1=Rt((j-1)*n1+1:(j-1)*n1+n1);
        d(i,j)=A*Rt1;
        if mod(j,2)~=0
            E(k)=E(k)+d(i,j);
        else
            E(k)=E(k)+d(i,j);
        end
    end
    k=k+1;
    if k>L
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=T+G;%二阶差分统计量的平移量

%%%%%%%%%%%%%%%%%%嵌入水印%%%%%%%%%%%%%%%%%%%%%
for j=1:S/n1
    bz(j)=floor((B+(j-1))/(S/n1));  %计算二阶差分统计量为正数时嵌入水印过程中每个二阶差分的变化量
    bf(j)=floor((-B+(j-1))/(S/n1)); %计算二阶差分统计量为负数时嵌入水印过程中每个二阶差分的变化量
end
% for j=1:S/3
%     bz(j)=floor((B+(j-1))/(S/3))*(-1)^(mod(j,2)+1);  %计算二阶差分统计量为正数时嵌入水印过程中每个二阶差分的变化量
%     bf(j)=floor((-B+(j-1))/(S/3))*(-1)^(mod(j,2)+1); %计算二阶差分统计量为负数时嵌入水印过程中每个二阶差分的变化量
% end

n2=2^n; %每个n阶差分的变化量需要平分的个数

Kmax=ceil((floor((B+(S/n1-1))/(S/n1))-floor((floor((B+(S/n1-1))/(S/n1)))/n2)*(2^n-2))/2); %嵌入系数
Kmin=floor(floor(B/(S/n1))/n2);

k=1;
for i=start:L1
    Rt=Rw((i-1)*S+1:(i-1)*S+S);
    if wm(k)==1&&E(k)>=0&&E(k)<T
        for j=1:S/n1
          Rt1=Rt((j-1)*n1+1:(j-1)*n1+n1);
          b1=floor(bz(j)/n2);
          b3=0;
          for kk=2:n1-1
              Rt1(kk)=Rt1(kk)+b1*(-1)^(kk+1);
              b3=b3+A1(kk)*b1;
          end
          b2=bz(j)-b3;
          Rt1(1)=Rt1(1)+floor(b2/2);
          Rt1(n1)=Rt1(n1)+floor((b2+1)/2)*(-1)^(n1+1);
          Rt((j-1)*n1+1:(j-1)*n1+n1)=Rt1;
        end
    elseif wm(k)==1&&E(k)<0&&E(k)>-T
        for j=1:S/n1
          Rt1=Rt((j-1)*n1+1:(j-1)*n1+n1);
          b1=ceil(bf(j)/n2);
          b3=0;
          for kk=2:n1-1
              Rt1(kk)=Rt1(kk)+b1*(-1)^(kk+1);
              b3=b3+A1(kk)*b1;
          end
          b2=bf(j)-b3;
          Rt1(1)=Rt1(1)+floor(b2/2);
          Rt1(n1)=Rt1(n1)+floor((b2+1)/2)*(-1)^(n1+1);
          Rt((j-1)*n1+1:(j-1)*n1+n1)=Rt1;
        end
    end
    Rw((i-1)*S+1:(i-1)*S+S)=Rt;      
    k=k+1;
    if k>L
        break;
    end
end





