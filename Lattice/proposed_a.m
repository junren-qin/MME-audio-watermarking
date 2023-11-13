function a=proposed_a(B_name,R)
    [B,rp,rc,G]=lattice_information(B_name);
    N=size(B,1);
    V_sphere=((pi^(N/2))*(rp^N))/gamma(1+N/2);
    a=1-(V_sphere/det(R*B))^(1/N);
end