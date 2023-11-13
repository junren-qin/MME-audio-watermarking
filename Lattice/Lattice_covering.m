function m=Lattice_covering(a)
    switch a
    case '1'
        m=1;
    case 'Z'
        m=sqrt(2);
    case 'A2'
        m=2/sqrt(3);
    case 'A3'
        m=sqrt(2);
    case 'D4'
        m=sqrt(2);
    case 'E8'
        m=sqrt(2);
    end
end