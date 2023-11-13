function m=Lattice_Basis(a,N)
switch a
    case 'Z'
        if nargin==1
            dist("please input dimension of Z lattice");
        else
            m=eye(N);
        end
    case 'A2'
        m=[[sqrt(3)/2;.5],[0;1]];
    case 'A3-self'
        m=[[1;0;0],[1/2;sqrt(3)/2;0],[1/2;sqrt(3)/6;sqrt(6)/3]];
    case 'A3'
        m=[[-1;-1;0],[1;-1;0],[0;1;-1]];
    case 'A3*'
        m=[[2;0;0],[0;2;0],[1;1;1]];
    case 'D4'
%         m=[[eye(3);zeros(1,3)],ones(4,1)];
        m=eye(4)+[ones(1,4);zeros(3,4)];
    case 'D5'
        m=[[eye(4);zeros(1,4)],ones(5,1)/2];
    case 'E6'
        m=[[zeros(1,5);-1*eye(5);zeros(2,5)]+[zeros(2,5);eye(5);zeros(1,5)],[ones(4,1)*1/2;ones(4,1)*(-1/2)]];
    case 'E7'
        m=[[-1*eye(6);zeros(2,6)]+[zeros(1,6);eye(6);zeros(1,6)],[ones(4,1)*(1/2);ones(4,1)*(-1/2)]];
    case 'E8'
        m=[[2;zeros(7,1)],[-1*eye(6);zeros(2,6)]+[zeros(1,6);eye(6);zeros(1,6)],(1/2)*ones(8,1)];
end