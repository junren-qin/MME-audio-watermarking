%%%%%%%%%%%%%%%%%%��ˮӡ�ӳ���%%%%%%%%%%%%%%%%%%%%
function [ Rc,wr,Ew,dw ] = Detecting( Rw,S,T,G,L,n)

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
A1=abs(A);%����n�ײ�ַ��̵��޷���ϵ������ȡ���̹�����Ҫ�õ�

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
for k=1:L
    if (Ew(k)<T) && (Ew(k)>-T)
        wr(k)=0;
    else
        wr(k)=1;
    end
end

%%%%%%%%%%%%%%%%%%�ָ���Ƶ%%%%%%%%%%%%%%%%%%%%%
for j=1:S/n1
    bz(j)=floor((B+(j-1))/(S/n1));  %������ײ��ͳ����Ϊ����ʱǶ��ˮӡ������ÿ�����ײ�ֵı仯��
    bf(j)=floor((-B+(j-1))/(S/n1)); %������ײ��ͳ����Ϊ����ʱǶ��ˮӡ������ÿ�����ײ�ֵı仯��
end
% for j=1:S/3
%     bz(j)=floor((B+(j-1))/(S/3))*(-1)^(mod(j,2)+1);  %������ײ��ͳ����Ϊ����ʱǶ��ˮӡ������ÿ�����ײ�ֵı仯��
%     bf(j)=floor((-B+(j-1))/(S/3))*(-1)^(mod(j,2)+1); %������ײ��ͳ����Ϊ����ʱǶ��ˮӡ������ÿ�����ײ�ֵı仯��
% end

n2=2^n; %ÿ��n�ײ�ֵı仯����Ҫƽ�ֵĸ���

k=1;
for i=start:L1
    Rt=Rc((i-1)*S+1:(i-1)*S+S);
    if Ew(k)>=T+G
        for j=1:S/n1
          Rt1=Rt((j-1)*n1+1:(j-1)*n1+n1);
          b1=floor(bz(j)/n2);
          b3=0;
          for kk=2:n1-1
              Rt1(kk)=Rt1(kk)-b1*(-1)^(kk+1);
              b3=b3+A1(kk)*b1;
          end
          b2=bz(j)-b3;
          Rt1(1)=Rt1(1)-floor(b2/2);
          Rt1(n1)=Rt1(n1)-floor((b2+1)/2)*(-1)^(n1+1);
          Rt((j-1)*n1+1:(j-1)*n1+n1)=Rt1;
        end
    elseif Ew(k)<=-(T+G)
        for j=1:S/n1
          Rt1=Rt((j-1)*n1+1:(j-1)*n1+n1);
          b1=ceil(bf(j)/n2);
          b3=0;
          for kk=2:n1-1
              Rt1(kk)=Rt1(kk)-b1*(-1)^(kk+1);
              b3=b3+A1(kk)*b1;
          end
          b2=bf(j)-b3;
          Rt1(1)=Rt1(1)-floor(b2/2);
          Rt1(n1)=Rt1(n1)-floor((b2+1)/2)*(-1)^(n1+1);
          Rt((j-1)*n1+1:(j-1)*n1+n1)=Rt1;
        end
    end
    Rc((i-1)*S+1:(i-1)*S+S)=Rt;    
    k=k+1;
    if k>L
        break;
    end  
end



