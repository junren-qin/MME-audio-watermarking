%%%%%%%%%%%%%%%%%%��ˮӡ�ӳ���%%%%%%%%%%%%%%%%%%%%
function [ wr,wr1,wr2,wr3,Ew,dw ] = Detecting_A( Rw,S,T,G,L,n)

Rc=Rw;  
Lg=length(Rc);     %���������
L1=floor(Lg/S);     %L1��֡
n1=n+1;
start=1;     %Ƕ�����

%%%%%%%%%%%%%%%%%%������ײ��ͳ����%%%%%%%%%%%%%%%
A=zeros(1,n1);
for i=1:n1
    A(i)=(-1)^(i+1)*nchoosek(n,i-1); %����n�ײ�ַ��̵Ĵ�����ϵ��
end

Ew=zeros(1,L);
k=1;
for i=start:L1
    Rt=Rc((i-1)*S+1:(i-1)*S+S);
    for j=1:S/n1
        Rt1=Rt((j-1)*n1+1:(j-1)*n1+n1);
        dw(i,j)=A*Rt1;
        if mod(j,2)~=0
            Ew(k)=Ew(k)+dw(i,j);
        else
            Ew(k)=Ew(k)+dw(i,j);
        end
    end
    k=k+1;
    if k>L
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=T+G;%���ײ��ͳ������ƽ����

%%%%%%%%%%%%%%%%%%��ȡˮӡ%%%%%%%%%%%%%%%%%%%%%
wr=zeros(1,L);
wr1=zeros(1,L);
wr2=zeros(1,L);
wr3=zeros(1,L);

%%%%%�о�1%%%%%%%
Q=T+G/10;
for k=1:L
    if (Ew(k)<Q) && (Ew(k)>-Q)
        wr1(k)=0;
    else
        wr1(k)=1;
    end
end

%%%%%�о�2%%%%%%%
Q=T+G*0.5;
for k=1:L
    if (Ew(k)<Q) && (Ew(k)>-Q)
        wr2(k)=0;
    else
        wr2(k)=1;
    end
end

%%%%%�о�3%%%%%%%
% [g,C]=kmeans(Ew,3);
[g,C]=kmeans(Ew',3,'Distance','sqEuclidean','Start','sample','Replicates',20);       %kmeans����

if C(1)==max(C)&&C(2)==min(C)
    k=3;
elseif C(1)==min(C)&&C(2)==max(C)
    k=3;
elseif C(1)==max(C)&&C(3)==min(C)
    k=2;
elseif C(1)==min(C)&&C(3)==max(C)
    k=2;
elseif C(2)==max(C)&&C(3)==min(C)
    k=1;
elseif C(2)==min(C)&&C(3)==max(C)
    k=1;
end

for i=1:L
    if g(i)==k
        wr3(i)=0;      
    else 
        wr3(i)=1;        
    end
end

%%%%����3���о�����ȡˮӡ%%%%%%%
for k=1:L
    if (wr1(k)==wr2(k))&&(wr1(k)==wr3(k))
        wr(k)=wr1(k);
    elseif (wr1(k)==wr2(k))&&(wr1(k)~=wr3(k))
        wr(k)=wr1(k);
    elseif (wr2(k)==wr3(k))&&(wr1(k)~=wr2(k))
        wr(k)=wr2(k);
    elseif (wr1(k)==wr3(k))&&(wr1(k)~=wr2(k))
        wr(k)=wr1(k);
    end
end







