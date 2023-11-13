function r=Lattice_packing(a)
% switch a
%     case 'Z'
%         r=1/2;
%     case 'A2'
%         r=1/2;
%     case 'D4'
%         r=1/sqrt(2);
%     case 'E8'
%         r=1/sqrt(2);
%     B=Lattice_Basis(a);
    B=a;
    N=size(B,1);
    k=zeros(1,N);
    for i=1:N
        k(i)=norm(B(:,i),2);
    end
    r=min(k)/2;
end