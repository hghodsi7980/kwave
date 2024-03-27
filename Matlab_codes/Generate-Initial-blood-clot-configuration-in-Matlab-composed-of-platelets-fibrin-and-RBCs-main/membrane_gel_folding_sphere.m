function [ Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang,Tdih,N1dih,N2dih,N3dih,N4dih, output111, output11, output1, output1111, output2222, output222] = ...
    membrane_gel_folding_sphere( Smap, Box, initial_number_points, Leq, ave_connectivity, r_cutoff, H, ~, MatomOne, periodic_x, periodic_y, periodic_z, atom_type, mol_type, bond_type, angle_type, psize, gelsize, percentActive,...
    pdata, crosslinkData, p_cutoff,maxConnections,num_p,r_cutoff_inner,clot_sidelength,num_rbc,box_sidelength,data_coords_p,bonds_p)

Cx = [];
Cy = [];
Cz = [];
Natom = [];
Matom = [];
Tmol = [];
Tatom = [];
Tbond = [];
N1bond = [];
N2bond = [];
Tang = [];
N1ang = [];
N2ang = [];
N3ang = [];

% cutoff distance - periodic - more organized - new point, bonds and angles
L.x = gelsize;L.y = H;L.z = gelsize;
% L.x = gelsize;L.y = H;L.z = gelsize;

%initial_number_points = round(initial_number_points * L.x * L.y * L.z);
initial_number_points = length(crosslinkData);

output1 = [];           %initial points

% initial_number_points = round ( (L.x * L.y * L.z) * Initial_Network_Rho );

% output1 = textread(sprintf('nodes_%d_1.dat', initial_number_points));
% output1 = textread(sprintf('nodes_%d_1.dat', 1000));
% output1(:,2) = atom_type;
% output1(:,6) = mol_type;
output1(:,1) = crosslinkData(:,1);
output1(:,2) = crosslinkData(:,3);
output1(:,3) = crosslinkData(:,4);
output1(:,4) = crosslinkData(:,5);
output1(:,5) = crosslinkData(:,6);
output1(:,6) = crosslinkData(:,2);

output11 = [];
output111 = [];          %final points
output2 = [];           %initial bonds
output22 = [];
output222 = [];          %final bonds
output33 = [];          %angles

initial_connections = [];
global_bonds_counter = 0;
local_bonds_counter = zeros(length(output1(:,1)),1);
distance = 0;
number_divisions = 0;
added_points_counter_1 = 0;
added_points_counter_2 = 0;
new_bonds_counter_1 = 0;
new_bonds_counter_2 = 0;
angles_counter = 0;

rg_molecule = mol_type;

for i=1:length(output1(:,1))
    if ( local_bonds_counter(i,1) >= ave_connectivity )
        continue;
    end

    local_connections = zeros ( (ave_connectivity - local_bonds_counter(i,1)) , 2 );
    local_connections(:,1) = local_connections(:,1) + 2 * sqrt(L.x^2+L.y^2+L.z^2) ;

    for j=1:length(output1(:,1))
        if (j <= i)
            continue;
        end

        [C,I] = max ( local_connections(:,1) );

        distance = sqrt ( ( min( abs(output1(i,3) - output1(j,3)), L.x - abs(output1(i,3) - output1(j,3)) ) )^2 ...
            + ( min( abs(output1(i,4) - output1(j,4)), L.y - abs(output1(i,4) - output1(j,4)) ) )^2 ...
            + ( min( abs(output1(i,5) - output1(j,5)), L.z - abs(output1(i,5) - output1(j,5)) ) )^2 );

        if ( ( distance < C ) && ( local_bonds_counter(j,1) < ave_connectivity ) && ( distance < r_cutoff )&& ( distance > r_cutoff_inner))
            local_connections(I,1) = distance;
            local_connections(I,2) = j;
        end
    end
    initial_connections ( i , ( local_bonds_counter(i,1) + 3 ) : ave_connectivity + 2) = local_connections(:,2)';

    for k=1:length( local_connections(:,2) )
        if ( local_connections(k,2) == 0 )
            continue;
        end
        global_bonds_counter = global_bonds_counter + 1;
        local_bonds_counter(i,1) = local_bonds_counter(i,1) +1;
        local_bonds_counter( local_connections(k,2),1) = local_bonds_counter( local_connections(k,2),1 ) + 1;
        initial_connections ( local_connections(k,2) , local_bonds_counter( local_connections(k,2) ) + 2 ) = i;
        output2( global_bonds_counter, 1) = global_bonds_counter;
        output2( global_bonds_counter, 2) = bond_type;
        output2( global_bonds_counter, 3) = i;
        output2( global_bonds_counter, 4) = local_connections(k,2);
        output2( global_bonds_counter, 5) = local_connections(k,1);
    end

end

initial_connections(:,1) = (1:length(output1(:,1)))';
initial_connections(:,2) = local_bonds_counter(:,1);
avg_fib_length = mean(output2(:,end));


for i=1:length(output1(:,1))
    output11(i,:) = output1(i,:);
end

for ii=1:length(output2(:,1))
    flagx = 0;
    flagy = 0;
    flagz = 0;
    if ( abs( output1(output2(ii,3),3) - output1( output2(ii,4),3) ) > L.x/2 )
        flagx = 1;
    end
    if ( abs( output1(output2(ii,3),4) - output1( output2(ii,4),4) ) > L.y/2 )
        flagy = 1;
    end
    if ( abs( output1(output2(ii,3),5) - output1( output2(ii,4),5) ) > L.z/2 )
        flagz = 1;
    end


    if ( ( (flagx == 1) && (periodic_x==0) ) || ( (flagy == 1) && (periodic_y==0) ) || ( (flagz == 1) && (periodic_z==0) ) )

        added_points_counter_1 = added_points_counter_1 + 2;
        added_points = zeros(8,3);

        added_points(1,1) = output1( output2(ii,4),3 ) - sign(output1( output2(ii,4),3) ) * L.x * flagx;
        added_points(1,2) = output1( output2(ii,4),4 ) - sign(output1( output2(ii,4),4) ) * L.y * flagy;
        added_points(1,3) = output1( output2(ii,4),5 ) - sign(output1( output2(ii,4),5) ) * L.z * flagz;

        added_points(5,1) = output1( output2(ii,3),3 ) - sign(output1( output2(ii,3),3) ) * L.x * flagx;
        added_points(5,2) = output1( output2(ii,3),4 ) - sign(output1( output2(ii,3),4) ) * L.y * flagy;
        added_points(5,3) = output1( output2(ii,3),5 ) - sign(output1( output2(ii,3),5) ) * L.z * flagz;



        if ( (periodic_x == 0) && ( flagx == 1 ) )

            if ( abs(added_points(1,1)) > (L.x/2) )
                added_points(2,1) = sign( added_points(1,1) ) * L.x / 2;
            else
                added_points(2,1) = added_points(1,1);
            end

            added_points(2,2) = output1( output2(ii,3),4 ) + ...
                ( ( added_points(1,2) - output1( output2(ii,3),4 ) ) / ( added_points(1,1) - output1( output2(ii,3),3 ) ) ) ...
                * ( added_points(2,1) - output1( output2(ii,3),3 ) );
            added_points(2,3) = output1( output2(ii,3),5 ) + ...
                ( ( added_points(1,3) - output1( output2(ii,3),5 ) ) / ( added_points(1,1) - output1( output2(ii,3),3 ) ) ) ...
                * ( added_points(2,1) - output1( output2(ii,3),3 ) );

            if ( abs(added_points(5,1)) > (L.x/2) )
                added_points(6,1) = sign( added_points(5,1) ) * L.x / 2;
            else
                added_points(6,1) = added_points(5,1);
            end

            added_points(6,2) = output1( output2(ii,4),4 ) + ...
                ( ( added_points(5,2) - output1( output2(ii,4),4 ) ) / ( added_points(5,1) - output1( output2(ii,4),3 ) ) ) ...
                * ( added_points(6,1) - output1( output2(ii,4),3 ) );
            added_points(6,3) = output1( output2(ii,4),5 ) + ...
                ( ( added_points(5,3) - output1( output2(ii,4),5 ) ) / ( added_points(5,1) - output1( output2(ii,4),3 ) ) ) ...
                * ( added_points(6,1) - output1( output2(ii,4),3 ) );
        else
            added_points(2,1) = added_points(1,1);
            added_points(2,2) = added_points(1,2);
            added_points(2,3) = added_points(1,3);

            added_points(6,1) = added_points(5,1);
            added_points(6,2) = added_points(5,2);
            added_points(6,3) = added_points(5,3);
        end


        if ( (periodic_y == 0) && ( flagy == 1 ) )

            if ( abs(added_points(2,2)) > (L.y/2) )
                added_points(3,2) = sign( added_points(2,2) ) * L.y / 2;
            else
                added_points(3,2) = added_points(2,2);
            end

            added_points(3,1) = output1( output2(ii,3),3 ) + ...
                ( ( added_points(2,1) - output1( output2(ii,3),3 ) ) / ( added_points(2,2) - output1( output2(ii,3),4 ) ) ) ...
                * ( added_points(3,2) - output1( output2(ii,3),4 ) );
            added_points(3,3) = output1( output2(ii,3),5 ) + ...
                ( ( added_points(2,3) - output1( output2(ii,3),5 ) ) / ( added_points(2,2) - output1( output2(ii,3),4 ) ) ) ...
                * ( added_points(3,2) - output1( output2(ii,3),4 ) );

            if ( abs(added_points(6,2)) > (L.y/2) )
                added_points(7,2) = sign( added_points(6,2) ) * L.y / 2;
            else
                added_points(7,2) = added_points(6,2);
            end

            added_points(7,1) = output1( output2(ii,4),3 ) + ...
                ( ( added_points(6,1) - output1( output2(ii,4),3 ) ) / ( added_points(6,2) - output1( output2(ii,4),4 ) ) ) ...
                * ( added_points(7,2) - output1( output2(ii,4),4 ) );
            added_points(7,3) = output1( output2(ii,4),5 ) + ...
                ( ( added_points(6,3) - output1( output2(ii,4),5 ) ) / ( added_points(6,2) - output1( output2(ii,4),4 ) ) ) ...
                * ( added_points(7,2) - output1( output2(ii,4),4 ) );
        else
            added_points(3,1) = added_points(2,1);
            added_points(3,2) = added_points(2,2);
            added_points(3,3) = added_points(2,3);

            added_points(7,1) = added_points(6,1);
            added_points(7,2) = added_points(6,2);
            added_points(7,3) = added_points(6,3);
        end




        if ( (periodic_z == 0) && ( flagz == 1 ) )

            if ( abs(added_points(3,3)) > (L.z/2) )
                added_points(4,3) = sign( added_points(3,3) ) * L.z / 2;
            else
                added_points(4,3) = added_points(3,3);
            end

            added_points(4,1) = output1( output2(ii,3),3 ) + ...
                ( ( added_points(3,1) - output1( output2(ii,3),3 ) ) / ( added_points(3,3) - output1( output2(ii,3),5 ) ) ) ...
                * ( added_points(4,3) - output1( output2(ii,3),5 ) );
            added_points(4,2) = output1( output2(ii,3),4 ) + ...
                ( ( added_points(3,2) - output1( output2(ii,3),4 ) ) / ( added_points(3,3) - output1( output2(ii,3),5) ) ) ...
                * ( added_points(4,3) - output1( output2(ii,3),5 ) );

            if ( abs(added_points(7,3)) > (L.z/2) )
                added_points(8,3) = sign( added_points(7,3) ) * L.z / 2;
            else
                added_points(8,3) = added_points(7,3);
            end

            added_points(8,1) = output1( output2(ii,4),3 ) + ...
                ( ( added_points(7,1) - output1( output2(ii,4),3 ) ) / ( added_points(7,3) - output1( output2(ii,4),5 ) ) ) ...
                * ( added_points(8,3) - output1( output2(ii,4),5 ) );
            added_points(8,2) = output1( output2(ii,4),4 ) + ...
                ( ( added_points(7,2) - output1( output2(ii,4),4 ) ) / ( added_points(7,3) - output1( output2(ii,4),5 ) ) ) ...
                * ( added_points(8,3) - output1( output2(ii,4),5 ) );
        else
            added_points(4,1) = added_points(3,1);
            added_points(4,2) = added_points(3,2);
            added_points(4,3) = added_points(3,3);

            added_points(8,1) = added_points(7,1);
            added_points(8,2) = added_points(7,2);
            added_points(8,3) = added_points(7,3);
        end

        output11(length(output1(:,1)) + added_points_counter_1 - 1,1) = length(output1(:,1)) + added_points_counter_1 - 1;
        output11(length(output1(:,1)) + added_points_counter_1 - 1,2) = atom_type;
        output11(length(output1(:,1)) + added_points_counter_1 - 1,6) = mol_type;

        if ( abs(added_points(4,1)) > (L.x/2) )
            output11(length(output1(:,1)) + added_points_counter_1 - 1,3) = added_points(4,1) - sign(added_points(4,1)) * L.x;
        else
            output11(length(output1(:,1)) + added_points_counter_1 - 1,3) = added_points(4,1);
        end

        if ( abs(added_points(4,2)) > (L.y/2) )
            output11(length(output1(:,1)) + added_points_counter_1 - 1,4) = added_points(4,2) - sign(added_points(4,2)) * L.y;
        else
            output11(length(output1(:,1)) + added_points_counter_1 - 1,4) = added_points(4,2);
        end

        if ( abs(added_points(4,3)) > (L.z/2) )
            output11(length(output1(:,1)) + added_points_counter_1 - 1,5) = added_points(4,3) - sign(added_points(4,3)) * L.z;
        else
            output11(length(output1(:,1)) + added_points_counter_1 - 1,5) = added_points(4,3);
        end


        output11(length(output1(:,1)) + added_points_counter_1,1) = length(output1(:,1)) + added_points_counter_1;
        output11(length(output1(:,1)) + added_points_counter_1,2) = atom_type;
        output11(length(output1(:,1)) + added_points_counter_1,6) = mol_type;

        if ( abs(added_points(8,1)) > (L.x/2) )
            output11(length(output1(:,1)) + added_points_counter_1,3) = added_points(8,1) - sign(added_points(8,1)) * L.x;
        else
            output11(length(output1(:,1)) + added_points_counter_1,3) = added_points(8,1);
        end

        if ( abs(added_points(8,2)) > (L.y/2) )
            output11(length(output1(:,1)) + added_points_counter_1,4) = added_points(8,2) - sign(added_points(8,2)) * L.y;
        else
            output11(length(output1(:,1)) + added_points_counter_1,4) = added_points(8,2);
        end

        if ( abs(added_points(8,3)) > (L.z/2) )
            output11(length(output1(:,1)) + added_points_counter_1,5) = added_points(8,3) - sign(added_points(8,3)) * L.z;
        else
            output11(length(output1(:,1)) + added_points_counter_1,5) = added_points(8,3);
        end

        new_bonds_counter_1 = new_bonds_counter_1 + 2;
        output22(new_bonds_counter_1 - 1,1) = new_bonds_counter_1 - 1;
        output22(new_bonds_counter_1 - 1,2) = bond_type;
        output22(new_bonds_counter_1 - 1,3) = output2(ii,3);
        output22(new_bonds_counter_1 - 1,4) = length(output1(:,1)) + added_points_counter_1 - 1;
        output22(new_bonds_counter_1 - 1,5) = norm ( output1(output2(ii,3),3:5) - added_points(4,:) );

        output22(new_bonds_counter_1,1) = new_bonds_counter_1;
        output22(new_bonds_counter_1,2) = bond_type;
        output22(new_bonds_counter_1,3) = output2(ii,4);
        output22(new_bonds_counter_1,4) = length(output1(:,1)) + added_points_counter_1;
        output22(new_bonds_counter_1,5) = norm ( output1(output2(ii,4),3:5) - added_points(8,:) );

    else

        new_bonds_counter_1 = new_bonds_counter_1 + 1;
        output22(new_bonds_counter_1,1) = new_bonds_counter_1;
        output22(new_bonds_counter_1,2) = bond_type;
        output22(new_bonds_counter_1,3) = output2(ii,3);
        output22(new_bonds_counter_1,4) = output2(ii,4);
        output22(new_bonds_counter_1,5) = output2(ii,5);
    end

end


for i=1:length(output11(:,1))
    output111(i,:) = output11(i,:);
end

vecArr = [];

for i=1:length(output22(:,1))
    number_divisions = round( output22(i,5) / Leq );
    vecArr = [vecArr; number_divisions];
    if ( number_divisions <= 1)
        new_bonds_counter_2 = new_bonds_counter_2 + 1;
        output222(new_bonds_counter_2,1) = new_bonds_counter_2;
        output222(new_bonds_counter_2,2) = bond_type;
        output222(new_bonds_counter_2,3) = output22(i,3);
        output222(new_bonds_counter_2,4) = output22(i,4);
        output222(new_bonds_counter_2,5) = output22(i,5);
        continue;
    end

    for j=1:number_divisions

        %New ponits
        if (j < number_divisions)
            added_points_counter_2 = added_points_counter_2 +1;

            output111(length(output11(:,1)) + added_points_counter_2,1) = length(output11(:,1)) + added_points_counter_2;
            output111(length(output11(:,1)) + added_points_counter_2,2) = atom_type;

            if ( abs( output11(output22(i,3),3) - output11(output22(i,4),3) ) > (L.x/2) )

                xdis = output11(output22(i,3),3) ...
                    + j * sign(output11(output22(i,3),3)) * ( L.x - abs( output11(output22(i,3),3) - output11(output22(i,4),3) ) ) / number_divisions;

                if ( xdis < L.x/2 && xdis > -L.x/2)
                    output111(length(output11(:,1)) + added_points_counter_2,3) = xdis;
                else
                    output111(length(output11(:,1)) + added_points_counter_2,3) = xdis - sign(xdis) * L.x;
                end
            else
                output111(length(output11(:,1)) + added_points_counter_2,3) = output11(output22(i,3),3) ...
                    + j * ( output11(output22(i,4),3) - output11(output22(i,3),3) ) / number_divisions;
            end


            if ( abs( output11(output22(i,3),4) - output11(output22(i,4),4) ) > (L.y/2) )

                ydis = output11(output22(i,3),4)...
                    + j * sign(output11(output22(i,3),4)) *( L.y - abs ( output11(output22(i,3),4) - output11(output22(i,4),4) ) ) / number_divisions;

                if ( ydis < L.y/2 && ydis > -L.y/2)
                    output111(length(output11(:,1)) + added_points_counter_2,4) = ydis;
                else
                    output111(length(output11(:,1)) + added_points_counter_2,4) = ydis - sign(ydis) * L.y;
                end
            else
                output111(length(output11(:,1)) + added_points_counter_2,4) = output11(output22(i,3),4) ...
                    + j * ( output11(output22(i,4),4) - output11(output22(i,3),4) ) / number_divisions;
            end


            if ( abs( output11(output22(i,3),5) - output11(output22(i,4),5) ) > (L.z/2) )

                zdis = output11(output22(i,3),5)...
                    + j * sign(output11(output22(i,3),5)) * ( L.z - abs ( output11(output22(i,3),5) - output11(output22(i,4),5) ) ) / number_divisions;

                if ( zdis < L.z/2 && zdis > -L.z/2)
                    output111(length(output11(:,1)) + added_points_counter_2,5) = zdis;
                else
                    output111(length(output11(:,1)) + added_points_counter_2,5) = zdis - sign(zdis) * L.z;
                end
            else
                output111(length(output11(:,1)) + added_points_counter_2,5) = output11(output22(i,3),5) ...
                    + j * ( output11(output22(i,4),5) - output11(output22(i,3),5) ) / number_divisions;
            end

            if j==1
                rg_molecule = rg_molecule+1;
            end
        end

        if number_divisions ==1;
            output111(length(output11(:,1)) + added_points_counter_2,6) = mol_type;
        else
            output111(length(output11(:,1)) + added_points_counter_2,6) = rg_molecule;
        end


        %New bonds
        if (j == 1)
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = output22(i,3);
            output222(new_bonds_counter_2,4) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), L.x - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), L.y - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), L.z - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
        elseif (j == number_divisions)
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,4) = output22(i,4);
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), L.x - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), L.y - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), L.z - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
        else
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = length( output11(:,1) ) + added_points_counter_2 - 1;
            output222(new_bonds_counter_2,4) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), L.x - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), L.y - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), L.z - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
        end

        %New angles
        if ( (j == 1) && (number_divisions == 2) )
            angles_counter = angles_counter + 1;
            output33(angles_counter,1) = angles_counter;
            output33(angles_counter,2) = angle_type;
            output33(angles_counter,3) = output22(i,3);
            output33(angles_counter,4) = length( output11(:,1) ) + added_points_counter_2;
            output33(angles_counter,5) = output22(i,4);
        elseif ( (j == 1) && (number_divisions ~= 2) )
            angles_counter = angles_counter + 1;
            output33(angles_counter,1) = angles_counter;
            output33(angles_counter,2) = angle_type;
            output33(angles_counter,3) = output22(i,3);
            output33(angles_counter,4) = length( output11(:,1) ) + added_points_counter_2;
            output33(angles_counter,5) = length( output11(:,1) ) + added_points_counter_2 + 1;
        elseif ( (j == number_divisions - 1) && (number_divisions ~= 2) )
            angles_counter = angles_counter + 1;
            output33(angles_counter,1) = angles_counter;
            output33(angles_counter,2) = angle_type;
            output33(angles_counter,3) = length( output11(:,1) ) + added_points_counter_2 - 1;
            output33(angles_counter,4) = length( output11(:,1) ) + added_points_counter_2;
            output33(angles_counter,5) = output22(i,4);
        elseif ( (1 < j) && ( j < number_divisions - 1) )
            angles_counter = angles_counter + 1;
            output33(angles_counter,1) = angles_counter;
            output33(angles_counter,2) = angle_type;
            output33(angles_counter,3) = length( output11(:,1) ) + added_points_counter_2 - 1;
            output33(angles_counter,4) = length( output11(:,1) ) + added_points_counter_2;
            output33(angles_counter,5) = length( output11(:,1) ) + added_points_counter_2 + 1;

        end

    end

end

atom_flag = zeros ( length(output111(:,1)),2 );
bond_flag = zeros ( length(output222(:,1)),1 );
angle_flag = zeros ( length(output33(:,1)),1 );
atom_flag_counter = 0;
bond_flag_counter = 0;
angle_flag_counter = 0;

output1111 = [];
output11111 = []; % I added this for DFS
output2222 = [];
output22222 = [];
output33333 = [];
output333 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BELOW FOR DFS TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(output111(:,1))
    %     if ( sqrt( output111(i,3)^2 + output111(i,5)^2 + output111(i,4)^2) - (L.x./2) >= 0)
    %         continue;
    %     end
    atom_flag_counter = atom_flag_counter + 1;
    atom_flag(i,1) = 1;
    atom_flag(i,2) = atom_flag_counter;
    output11111(atom_flag_counter,1) = atom_flag_counter;
    output11111(atom_flag_counter,3) = output111(i,3);
    output11111(atom_flag_counter,4) = output111(i,4);
    output11111(atom_flag_counter,5) = output111(i,5);
    output11111(atom_flag_counter,6) = mol_type;
    output11111(atom_flag_counter,2) = atom_type;
end

for i=1:length(output222(:,1))
    if ( (atom_flag(output222(i,3),1) == 1) && ( atom_flag(output222(i,4),1) == 1 ) )
        bond_flag(i) = 1;
        bond_flag_counter = bond_flag_counter + 1;
        output22222(bond_flag_counter,1) = bond_flag_counter;
        output22222(bond_flag_counter,2) = output222(i,2);
        output22222(bond_flag_counter,3) = atom_flag(output222(i,3),2);
        output22222(bond_flag_counter,4) = atom_flag(output222(i,4),2);
        output22222(bond_flag_counter,5) = output222(i,5);
        %         if (sqrt(output111(output22222(bond_flag_counter,3),3))^2 + output111(output22222(bond_flag_counter,3),4)^2 + output111(output22222(bond_flag_counter,3),5)^2) <= (gelsize./3);
        %             output22222(bond_flag_counter,2) = output22222(bond_flag_counter,2) + 1;
        %         end
    end
end

for i=1:length(output33(:,1))
    if ( (atom_flag(output33(i,3),1) == 1) && ( atom_flag(output33(i,4),1) == 1 ) && ( atom_flag(output33(i,5),1) == 1 ) )
        angle_flag(i) = 1;
        angle_flag_counter = angle_flag_counter + 1;
        output33333(angle_flag_counter,1) = angle_flag_counter;
        output33333(angle_flag_counter,2) = output33(i,2);
        output33333(angle_flag_counter,3) = atom_flag(output33(i,3),2);
        output33333(angle_flag_counter,4) = atom_flag(output33(i,4),2);
        output33333(angle_flag_counter,5) = atom_flag(output33(i,5),2);
    end
end

output1111 = [output111; [(length(output111) + (1:length(pdata))'), pdata(:,2:5), pdata(:,1)]];
avg_fibrinL = mean(output22(:,5));sprintf('average fibrin length=%.2fdpd',avg_fibrinL)
num_fibrin=length(output22);


%% add in platelets

num_per_p=length(data_coords_p);      bonds=bonds_p;
numentry = num_p *length(bonds);
BondinPlatelet = [];
for i = 1:num_p
    add = bonds;
    add(:,3)= bonds(:,3)+(i-1).*num_per_p;
    add(:,4)= bonds(:,4)+(i-1).*num_per_p;
    BondinPlatelet = [BondinPlatelet; add];
end
BondinPlatelet(:,1)=1:numentry;
BondinPlatelet(:,3) = BondinPlatelet(:,3)+ atom_flag_counter;
BondinPlatelet(:,4) = BondinPlatelet(:,4)+ atom_flag_counter;
output22222 = [output22222; [(length(output22222)+ (1:length(BondinPlatelet))'),BondinPlatelet(:,2:5)]];

save1 = output1111;
save2 = output22222;
save3 = output33333;

[my2,my1]=prepare20220113(output1111,output22222,H,gelsize,p_cutoff,maxConnections,num_per_p,num_p);
output1111 = my1;
output22222 = [output22222;my2];

% make fibrin bonds type 1
kkk=find(output22222(:,2)==2);
output22222(kkk,2)=1;

%% add plt initial filopodia
% final filopodia length = 12dpd (6 micron)
% initial filopodia length = 3 =1.5 in dpd length

N_seg=3;  % pltcenter--b1--b2--b3
L_seg=2;
pcid = find(output1111(:,2)>=4 & output1111(:,2)<=23);

% 10 groups:
% 9 groups:
group1type = (4:1:15)';  group2type = (16:1:23)'; group2type=[0;group2type];

filo_end_type = zeros(length(group2type), 3);
filo_end_type(1,:)=[24 24 25];
for i=1:length(filo_end_type(:,1))
    filo_end_type(i,:)=[24 24 25]+ (i-1).*[0 0 1];
end
filo_bond_types = filo_end_type;
for i=1:length(filo_bond_types(:,1))
    filo_bond_types(i,:)=[6 7 8]+ (i-1).*[1 1 1];
end

for i = 1:length(pcid)
    pc_type = output1111(pcid(i),2);
    if pc_type<=15
        group_id = 1;
    elseif pc_type >15
        group_id = find(group2type==pc_type);
    end
    myends = filo_end_type(group_id,3) .* ones(12,1) +[0; 0; 0; 0; 1; 1; 1;1; 2; 2; 2;2 ];
    bbb = filo_bond_types(group_id,1).* ones(12,1) +[0; 0; 0; 0; 1; 1; 1;1; 2; 2; 2 ;2];

    pcxyz = output1111(pcid(i),3:5);
    TH = 2*pi*rand(1,12);PH = asin(-1+2*rand(1,12)); [X,Y,Z] = sph2cart(TH,PH,1);
    % plot3(X,Y,Z,'*r') ; hold on; vv = [0,0,0];plot3(vv(1),vv(2),vv(3),'*k');axis equal
    add_atom = [];                  add_bond = [];                   add_ang = [];
    atomf = max(output1111(:,1));   bondf = max(output22222(:,1));   angf =  max(output33333(:,1));
    for j = 1:length(X)
        %   add atoms
        add = zeros(3,6);
        add(1,3:5)=pcxyz +[X(j) Y(j) Z(j)].*L_seg./3;
        add(2,3:5)=pcxyz +[X(j) Y(j) Z(j)].*L_seg.*2./3;
        add(3,3:5)=pcxyz +[X(j) Y(j) Z(j)].*L_seg;
        add(:,1)= atomf + (1:3)';
        add(:,2)=[24;24;myends(j)];
        add(:,6)=[24;24;myends(j)];
        add_atom = [add_atom;add];
        atomf = atomf+N_seg;
        % add bonds
        addb = zeros(3,5);
        addb(:,2)=bbb(j);
        addb(:,3)=[pcid(i); add(1,1);add(2,1) ];
        addb(:,4)=[add(1,1);add(2,1); add(3,1)];
        add_bond = [add_bond;addb];
        % add angles
        adda = zeros(2,5);   adda(:,2)=3;
        adda(:,3)=[pcid(i);  add(1,1)];   adda(:,4)=[add(1,1);add(2,1)];  adda(:,5)=[add(2,1);add(3,1)];
        add_ang = [add_ang;adda];
    end
    add_bond(:,1)= bondf+(1:length(add_bond(:,1)))';   add_ang(:,1)= angf+(1:length(add_ang(:,1)))';
    output1111=[output1111;add_atom];
    output22222=[output22222;add_bond];
    output33333=[output33333;add_ang];
end
sprintf('Total bead number fib and plt = %d',length(output1111(:,1)))

%% rbc

if num_rbc>0
    load('rbc_atomdata')
    load('rbc_bonddata')
    load('rbc_dihdata')
    save_atomdata_rbc=rbc_atomdata;  save_bonddata_rbc=rbc_bonddata;    save_dihedraldata_rbc= rbc_dihdata;

    % atom data
    rbc_shell_atomtype = 37;
    save_atomdata_rbc(:,2)=rbc_shell_atomtype; save_atomdata_rbc(:,3)=rbc_shell_atomtype;
    save_atomdata_rbc(:,4)= save_atomdata_rbc(:,4)- mean(save_atomdata_rbc(:,4));
    rbc_coords=save_atomdata_rbc;
    %     figure
    %     scatter3(rbc_coords(:,4),rbc_coords(:,5),rbc_coords(:,6),'b'); axis equal
    %     axis equal
    ss=0.28;
    rbc_coords(:,4) = rbc_coords(:,4)*ss;    rbc_coords(:,5) = rbc_coords(:,5)*ss;    rbc_coords(:,6) = rbc_coords(:,6)*ss;
    D_rbc_adj1 = max(rbc_coords(:,4))-min(rbc_coords(:,4));
    D_rbc_adj2 = max(rbc_coords(:,5))-min(rbc_coords(:,5));
    D_rbc_adj3 = max(rbc_coords(:,6))-min(rbc_coords(:,6));
    D_rbc_adj  = max([D_rbc_adj1;D_rbc_adj2;D_rbc_adj3]);
    bds=[];
    for i=1:length(save_bonddata_rbc(:,1))
        lsq = (rbc_coords(save_bonddata_rbc(i,3),4) - rbc_coords(save_bonddata_rbc(i,4),4) )^2 + ...
            (rbc_coords(save_bonddata_rbc(i,3),5) - rbc_coords(save_bonddata_rbc(i,4),5) )^2 +...
            (rbc_coords(save_bonddata_rbc(i,3),6) - rbc_coords(save_bonddata_rbc(i,4),6) )^2;
        bds(i)=sqrt(lsq);
    end
    bds=round(bds,3); bds=bds';
    sprintf('rbc shell atom type %d, rbc diameter=%.2fdpd, min bond L=%.4f, max bobd L=%.4f',rbc_shell_atomtype, D_rbc_adj,min(bds),max(bds) )


    %% generate rbc xyz  based on lattice, and then add in randomness
    xup = clot_sidelength /2+1; xlow = - clot_sidelength /2-1;
    yup = clot_sidelength /2+1; ylow = - clot_sidelength /2-1;
    zup = clot_sidelength /2+1; zlow = - clot_sidelength /2-1;


    nperside = ceil((num_rbc)^(1/3));   dprbc = (clot_sidelength+2)/nperside -0.1;
    xlist = (xlow+dprbc/2):dprbc:(xup);  xlist=xlist';    ylist = (ylow+dprbc/2):dprbc:(yup);  ylist=ylist';    zlist = (zlow+dprbc/2):dprbc:(zup);  zlist=zlist';
    tot_grid = length(xlist)*length(ylist)*length(zlist);
    rbc_xyz=[]; add=zeros(length(xlist),3);
    for qq =1:length(zlist)
        for uu = 1:length(ylist)
            add(:,1)=xlist;   add(:,2)=ylist(uu);   add(:,3)=zlist(qq);
            rbc_xyz=[rbc_xyz;add];
        end
    end
    ax = max(rbc_xyz(:,1)); bx = min(rbc_xyz(:,1));  mx = (ax+bx)/2;
    ay = max(rbc_xyz(:,2)); by = min(rbc_xyz(:,2));  my = (ay+by)/2;
    az = max(rbc_xyz(:,3)); bz = min(rbc_xyz(:,3));  mz = (az+bz)/2;
    rbc_xyz(:,1) = rbc_xyz(:,1)-mx;  rbc_xyz(:,2) = rbc_xyz(:,2)-my; rbc_xyz(:,3) = rbc_xyz(:,3)-mz;
    %   add in some randomness
    add_rand = -2 + 4.* rand(length(rbc_xyz),3);   rbc_xyz = rbc_xyz + add_rand;
    % make sure RBC center xyz do not overlay with plt and fibrin
    for bb = 1:length(rbc_xyz(:,1))
        cdistsq = (rbc_xyz(bb,1)-output1111(:,3)).^2+(rbc_xyz(bb,2)-output1111(:,4) ).^2+(rbc_xyz(bb,3)-output1111(:,5) ).^2;
        cdist = sqrt(cdistsq);   flagg = find(cdist<=(D_rbc_adj/2+0.3));
        iter=0; conti=1;
        if length(flagg)>=2
            while conti==1
                rbc_xyz(bb,:) =  rbc_xyz(bb,:) -1.5 + 3.* rand(1,3);         iter=iter+1;
                cdistsq2 = (rbc_xyz(bb,1)-output1111(:,3)).^2+(rbc_xyz(bb,2)-output1111(:,4) ).^2+(rbc_xyz(bb,3)-output1111(:,5) ).^2;
                cdist2 = sqrt(cdistsq2);   flagg2 = find(cdist2<=(D_rbc_adj/2+0.4));
                if length(flagg2)>=2
                    if iter>=4
                        conti=0;    rbc_xyz(bb,:)=0;
                    end
                end
                if length(flagg2)<2
                    conti=0;
                end
            end
        end
    end
    todelete = find(rbc_xyz(:,1)==0 & rbc_xyz(:,2)==0 & rbc_xyz(:,3)==0 );   N_del1=length(todelete);
    percdel = N_del1/num_rbc; 
    if percdel>=0.07
        sprintf('Warning! deleted %d out of %d RBC in clot', N_del1, num_rbc)
    else
        sprintf('deleted %d out of %d RBC xyz in clot', N_del1, num_rbc)
    end
    rbc_xyz(todelete,:)=[];
    final_num_rbc =length(rbc_xyz(:,1));
    while final_num_rbc>num_rbc
        delid = ceil(  ( length(rbc_xyz(:,1))-1) *rand(1)  );
        rbc_xyz(delid,:)=[];
        final_num_rbc =length(rbc_xyz(:,1));
    end


    %%
    RBC_atomdata=[];
    mol_ind = max(output1111(:,6));
    rbc_mol_start=mol_ind+1;
    sprintf('desired Nrbc is %d, actual N rbc is %d, rbc seperate by %.2fdpd, rbc mol start from %d', num_rbc, final_num_rbc,dprbc,rbc_mol_start )
    num_rbc=final_num_rbc;


    for i=1:num_rbc
        add_rbc = zeros(length(rbc_coords),5);
        add_rbc(:,2)=rbc_shell_atomtype;
        add_rbc(:,1)=rbc_mol_start;
        theta = rand.*360;   rbc_coordsr(:,4:6) = rbc_coords(:,4:6)*[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
        theta = rand.*360;   rbc_coordsrr(:,4:6) = rbc_coordsr(:,4:6)*[1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)];
        add_rbc(:,3)= rbc_xyz(i,1)+rbc_coordsrr(:,4);
        add_rbc(:,4)= rbc_xyz(i,2)+rbc_coordsrr(:,5);
        add_rbc(:,5)= rbc_xyz(i,3)+rbc_coordsrr(:,6);
        RBC_atomdata=[RBC_atomdata;add_rbc];
        rbc_mol_start=rbc_mol_start+1;
    end

    hRBC_atomdata=zeros(length(RBC_atomdata),6);
    hRBC_atomdata(:,1)=1:length(RBC_atomdata);
    hRBC_atomdata(:,2:6)=RBC_atomdata;
    output1111 = [output1111; [(length(output1111) + (1:length(hRBC_atomdata))'), hRBC_atomdata(:,3:6), hRBC_atomdata(:,2)]];

    % bond data
    num_bd_type = length(unique(save_bonddata_rbc(:,2)));
    ori_btype = unique(save_bonddata_rbc(:,2));
    save_bonddata_rbc(:,2)=save_bonddata_rbc(:,2)+14;
    cur_btype = unique(save_bonddata_rbc(:,2)); sprintf('rbc bond types are %d to %d',cur_btype(1),cur_btype(end))


    num_atom_per_rbc=length(rbc_coords);
    num_bond = num_rbc *length(save_bonddata_rbc);
    BondinRBC = [];
    for i = 1:num_rbc
        add = save_bonddata_rbc;
        add(:,3)= save_bonddata_rbc(:,3)+(i-1).*num_atom_per_rbc;
        add(:,4)= save_bonddata_rbc(:,4)+(i-1).*num_atom_per_rbc;
        BondinRBC = [BondinRBC; add];
    end
    BondinRBC(:,1)=1:num_bond;

    xxr = find(output1111(:,2)==rbc_shell_atomtype);
    atomflag = xxr(1)-1
    BondinRBC(:,3) = BondinRBC(:,3)+ atomflag;
    BondinRBC(:,4) = BondinRBC(:,4)+ atomflag;
    BondinRBC(:,5)=0.1;
    output22222 = [output22222; [(length(output22222)+ (1:length(BondinRBC))'),BondinRBC(:,2:5)]];

    % dihedral data
    rbc_di_type = unique(save_dihedraldata_rbc(:,2));
    rbc_dih = [];
    for i = 1:num_rbc
        add = save_dihedraldata_rbc;
        add(:,3)=save_dihedraldata_rbc(:,3)+(i-1).*num_atom_per_rbc;
        add(:,4)=save_dihedraldata_rbc(:,4)+(i-1).*num_atom_per_rbc;
        add(:,5)=save_dihedraldata_rbc(:,5)+(i-1).*num_atom_per_rbc;
        add(:,6)=save_dihedraldata_rbc(:,6)+(i-1).*num_atom_per_rbc;
        rbc_dih=[rbc_dih;add];
    end

    rbc_dih(:,3:6)=rbc_dih(:,3:6)+ atomflag;
    rbc_dih(:,1)=(1:length(rbc_dih));
    save_all_dihedral =  rbc_dih ;

end

%% add rbc innerfluid beads
type_inn = 38;
inner_numbead = 1100;
add_inner = [];
for i= 1:length(rbc_xyz(:,1))
    pos = rbc_xyz(i,:);
    add = zeros(inner_numbead,6);
    add(:,3)= pos(1)-0.1 +(0.2).*rand(inner_numbead,1);
    add(:,4)= pos(2)-0.1 +(0.2).*rand(inner_numbead,1);
    add(:,5)= pos(3)-0.1 +(0.2).*rand(inner_numbead,1);

    add_inner=[add_inner;add];
end
add_inner(:,1)= length(output1111(:,1))+ (1:length(add_inner(:,1)))';
add_inner(:,2)=type_inn;
add_inner(:,6)=type_inn;
output1111=[output1111;add_inner];


%% wall

% cap_len,cap_dia,num_rbc,clot_len)
%     xup = clot_sidelength /2; xlow = - clot_sidelength /2;
wall_dat=[];
addwall_dat=zeros(1,6);
addwall_dat(:,2)=39;  addwall_dat(:,6)=39;
addwall_dat(:,3)=xup+1; addwall_dat(:,4)=0; addwall_dat(:,5)=0;

wall_dat=[wall_dat;addwall_dat];

wall_dat2=[];
addwal2_dat=zeros(1,6);
addwal2_dat(:,2)=39;  addwal2_dat(:,6)=39;
addwal2_dat(:,3)=xlow-1; addwal2_dat(:,4)=0; addwal2_dat(:,5)=0;
wall_dat=[wall_dat;addwal2_dat];


output1111=[output1111;wall_dat];          output1111(:,1)=1:length(output1111(:,1));


%% change some fibrin type for different fix condition in flow

% allfib = find(output1111(:,2)==2); idmax = length(allfib)-1;
% select = sqrt(output1111(1:idmax,3).^2+output1111(1:idmax,4).^2 );
%
% attachmode = 2  % (1=only attach bottom; 2= all circular sides fixed)
% if attachmode==1
%     selected = find(select>(cap_dia/2-5) & output1111(1:idmax,3)<(min(output1111(1:idmax,3))+20) &  output1111(1:idmax,5)>-40 &  output1111(1:idmax,5)<40 );
% elseif attachmode==2
%     selected = find(select>(cap_dia/2-8) );
% end
% output1111(selected,2)=1;

%%
% figure
% k1 = find(output1111(:,2)==2); % fibrin
% figure
% scatter3(output1111(k1,3), output1111(k1,4), output1111(k1,5), '*b'   )
% hold on
% k1 = find(output1111(:,2)==1); % fibrin boundary
% scatter3(output1111(k1,3), output1111(k1,4), output1111(k1,5), '*k'   )
% axis equal
% hold on
%
% k2 = find(output1111(:,2)>=3 & output1111(:,2)<=24 ); % plt body
% scatter3(output1111(k2,3), output1111(k2,4), output1111(k2,5), '*g'   )
% hold on
%
% k2 = find(output1111(:,2)>=24 & output1111(:,2)<=34 ); % filopod
% scatter3(output1111(k2,3), output1111(k2,4), output1111(k2,5), '*g'   )
% hold on
% axis equal
%
% k4 = find(output1111(:,2)==35  );
% scatter3(output1111(k4,3), output1111(k4,4), output1111(k4,5), '*r'   )
% hold on
% axis equal
% k4 = find(output1111(:,2)==36  );
% scatter3(output1111(k4,3), output1111(k4,4), output1111(k4,5), '*r'   )
% hold on
% axis equal
%
%
% k4 = find(output1111(:,2)==38  ); % rbc
% scatter3(output1111(k4,3), output1111(k4,4), output1111(k4,5), '*k'   )
% hold on
% axis equal
% k4 = find(output1111(:,2)==37  ); % rbc
% scatter3(output1111(k4,3), output1111(k4,4), output1111(k4,5), '*r'   )
% hold on
% % k4 = find(output1111(:,2)==39  ); % vessel
% % scatter3(output1111(k4,3), output1111(k4,4), output1111(k4,5), '*y'   )
% % hold on
% axis equal


%%

Cx = output1111(:,3);
Cy = output1111(:,4);
Cz = output1111(:,5);
Natom = output1111(:,1);
Matom (1:length(output1111(:,1)))= MatomOne;
Tmol = output1111(:,6);
Tatom = output1111(:,2);


Tbond = output22222(:,2);
N1bond = output22222(:,3);
N2bond = output22222(:,4);


Tang = output33333(:,2);
N1ang = output33333(:,3);
N2ang = output33333(:,4);
N3ang = output33333(:,5);


Tdih =  save_all_dihedral(:,2);
N1dih = save_all_dihedral(:,3);
N2dih = save_all_dihedral(:,4);
N3dih = save_all_dihedral(:,5);
N4dih = save_all_dihedral(:,6);


end




