function create_output_file(FileName, Box, Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang,Tdih,N1dih,N2dih,N3dih,N4dih,box_sidelength)

Natom_tot = size(Natom);
Nbond_tot = size(Tbond);
Nang_tot = size(Tang);
Ndih_tot = size(Tdih);
Natom_tot = Natom_tot(1,1);
Nbond_tot = Nbond_tot(1,1);
Nang_tot = Nang_tot(1,1);
Ndih_tot = Ndih_tot(1,1);

if(Natom_tot > 0)
    Tatom_tot = max(Tatom)+1;  % for SRP bead
else
    Tatom_tot = 0;
end

if(Nbond_tot > 0)
    Tbond_tot = max(Tbond);
else
    Tbond_tot = 0;
end
 Tbond_tot=Tbond_tot+1
 
if(Nang_tot > 0)
    Tang_tot = max(Tang);
else
    Tang_tot = 0;
end

disp(sprintf('Total: Natom=%d  Nbond=%d  Nang=%d Tatom=%d  Tbond=%d  Tang=%d', ...
        Natom_tot, Nbond_tot, Nang_tot, Tatom_tot, Tbond_tot, Tang_tot));
disp(sprintf('Output file: %s', FileName));


fOutput = fopen(FileName, 'w+');

fprintf(fOutput, 'LAMMPS Description\n');
fprintf(fOutput, '\n');

fprintf(fOutput, '%d atoms\n', Natom_tot);
fprintf(fOutput, '%d bonds\n', Nbond_tot);
fprintf(fOutput, '%d angles\n', Nang_tot);
fprintf(fOutput, '%d dihedrals\n', Ndih_tot);
fprintf(fOutput, '\n');

fprintf(fOutput, '%d atom types\n', Tatom_tot);
fprintf(fOutput, '%d bond types\n', Tbond_tot);
fprintf(fOutput, '%d angle types\n', Tang_tot);
fprintf(fOutput, '363 dihedral types\n');

fprintf(fOutput, '\n');


xl=-box_sidelength/2; xh=box_sidelength/2; 
yl=-box_sidelength/2; yh=box_sidelength/2; 
zl=-box_sidelength/2; zh=box_sidelength/2; 

fprintf(fOutput, '%g %g xlo xhi\n', xl, xh);
fprintf(fOutput, '%g %g ylo yhi\n', yl, yh);
fprintf(fOutput, '%g %g zlo zhi\n', zl, zh);
fprintf(fOutput, '\n');

fprintf(fOutput, 'Masses\n');
fprintf(fOutput, '\n');

MatomT = zeros(Tatom_tot, 1);
for I=1:Natom_tot
    MatomT(Tatom(I)) = Matom(I);
end

for I=1:Tatom_tot
    if(MatomT(I) == 0)
        MatomT(I) = 1;
    end
    fprintf(fOutput, '%d %g\n', I, MatomT(I));
end
fprintf(fOutput, '\n');

fprintf(fOutput, 'Atoms\n');
fprintf(fOutput, '\n');
for I=1:Natom_tot
%    fprintf(fOutput, '%d %d %d %g %g %g\n', I, Tmol(I), Tatom(I),...
%    Cx(I), Cy(I), Cz(I)); %old version
% new version 07
    fprintf(fOutput, '%d %d %d %g %g %g\n', I, Tmol(I), Tatom(I), Cx(I), Cy(I), Cz(I));
end

fprintf(fOutput, '\n');

fprintf(fOutput, 'Bonds\n');
fprintf(fOutput, '\n');
for I=1:Nbond_tot
    fprintf(fOutput, '%d %d %d %d\n', I, Tbond(I), N1bond(I), N2bond(I));
end

fprintf(fOutput, '\n');

if(Nang_tot > 0)

    fprintf(fOutput, 'Angles\n');
    fprintf(fOutput, '\n');
    for I=1:Nang_tot
        fprintf(fOutput, '%d %d %d %d %d\n', I, Tang(I), N1ang(I), N2ang(I), N3ang(I));
    end
end

fprintf(fOutput, '\n');

if(Ndih_tot > 0)

    fprintf(fOutput, 'Dihedrals\n');
    fprintf(fOutput, '\n');
    for I=1:Ndih_tot
        fprintf(fOutput, '%d %d %d %d %d %d\n', I, Tdih(I), N1dih(I), N2dih(I), N3dih(I),N4dih(I));
    end
end



fclose(fOutput);

