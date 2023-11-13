function g=Lattice_Mean_Square(a)
    switch a
        case '1'
            g=1/12;
        case 'Z'
            g=1/12;
        case 'A2'
            g=5/(36*sqrt(3));
        case 'A3'
            g=0.078543;
        case 'D4'
            g=0.076603;
        case 'E8'
            g=0.071682;
    end
end