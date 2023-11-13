function [B,rp,rc,G]=lattice_information(lattice_name,N)
%     type=1350;
%     type=35*1e+3;
%     type=2000;
    type=2000;
    if lattice_name=='Z'
        if nargin==1
            N=2;
        end
        B=Lattice_Basis(lattice_name,N);
        rp=Lattice_packing(B);
        B=.5*B/rp*type;
        rp=.5*rp/rp*type;
        rc=sqrt(N)*rp*type;
    else
        B=Lattice_Basis(lattice_name);
        rp=Lattice_packing(B);
        B=.5*B/rp*type;
        rp=.5*rp/rp*type;
        rc=Lattice_covering(lattice_name)*rp*type;
    end
    G=Lattice_Mean_Square(lattice_name);   
end