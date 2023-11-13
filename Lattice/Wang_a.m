function a=Wang_a(B_name,R)
    [B,rp,rc,G]=lattice_information(B_name);
    a=(1-(rp/(R*rc)));
end