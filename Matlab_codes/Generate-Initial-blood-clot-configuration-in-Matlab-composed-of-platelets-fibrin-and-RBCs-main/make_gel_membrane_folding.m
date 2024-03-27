clear all
close all
clc

Lx =2500;
Ly = Lx;
Lz = Lx;

Box.x.lo = -Lx / 2.;
Box.x.hi = Lx / 2.;
Box.y.lo = -Ly / 2.;
Box.y.hi = Ly / 2.;
Box.z.lo = -Lz / 2.;
Box.z.hi = Lz / 2.;

L.x = Box.x.hi - Box.x.lo;L.y = Box.y.hi - Box.y.lo;L.z = Box.z.hi - Box.z.lo;

Rho = 3;
periodic_x = 0;periodic_y = 0;periodic_z = 0;
Nmap = 1;Smap = 2;
% Map = zeros(ceil(Smap*L.x+1), ceil(Smap*L.y+1), ceil(Smap*L.z+1));
MatomOne = 1;Mone = 1;
num_cluster = 200;R_particle = 0.209;
post_length = 40;post_inter_layer = 10;post_outer_layer = 11;post_angle = -0*pi/180;num_post = 1;
[Natom Matom Tmol Tatom Cx Cy Cz Tbond N1bond N2bond Tang N1ang N2ang N3ang] = create_arrays();

%%%%%%%%%%%%%% add gel_shell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gelsize = Lx/2;
gelsize = Lx./2;Initial_Network_Rho = 0.1;
percentActive = 0.5;Leq = 1;
p_cutoff=20;
maxConnections=0;

H = gelsize;
r1 = H - 1;
psize = 0.75;
atom_type = 2;
mol_type = 2;
bond_type = 2;
angle_type = 2;
%[Map NatomN MatomN TmolN TatomN CxN CyN CzN TbondN N1bondN N2bondN TangN N1angN N2angN N3angN] ...
%    = membrane_gel_folding_cross(Map, Smap, Box, Initial_Network_Rho, Leq, ave_connectivity, r_cutoff, H, r1, MatomOne, periodic_x, periodic_y, periodic_z, atom_type, mol_type, bond_type, angle_type, psize, gelsize);
% [Map NatomN MatomN TmolN TatomN CxN CyN CzN TbondN N1bondN N2bondN TangN N1angN N2angN N3angN] ...
%     = membrane_gel_folding_circle(Map, Smap, Box, Initial_Network_Rho, Leq, ave_connectivity, r_cutoff, H, r1, MatomOne, periodic_x, periodic_y, periodic_z, atom_type, mol_type, bond_type, angle_type, psize, gelsize);

% platelet data  (mol_type(3) atom_type(3) x y z)
% crosslink data (id mol_type(2) atom_type(2) x y z)

ave_connectivity = 5;
r_cutoff =44;
r_cutoff_inner = 15;

%% new crosslink and pdaata
clotsize_ref = 200;

num_cl_ref   = 5600;  % num_cl_ref= 1400;
num_p_ref    = 2000;
num_rbc_ref  = 3200;  % 

% set up ----   cubic box (rbc), sphere clot (fib plt)
clot_diameter = 50 ;  %(dpd)     
clot_sidelength = clot_diameter ; %(dpd)     50micron
box_sidelength = clot_sidelength+8;

sprintf('clot diameter=%ddpd, box size=%ddpd,c_cl=%d, c_plt=%d, c_rbc=%d',clot_diameter,box_sidelength,num_cl_ref,num_p_ref,num_rbc_ref)


cubicV = clot_sidelength^3 ;
boxV = box_sidelength^3 ;

clotV = (4/3)*pi()*(clot_diameter/2)^3;
volumeratio = clotV/cubicV;
num_rbc= ceil(num_rbc_ref *boxV/(clotsize_ref^3));

% calculating #fib crosslink, and generating location
num_CL= ceil(num_cl_ref *cubicV/(clotsize_ref^3));
xup = clot_sidelength /2; xlow =- clot_sidelength/2;   yup = clot_sidelength /2; ylow = -clot_sidelength/2;    zup = clot_sidelength/2;  zlow = -clot_sidelength/2;
new_crosslinkdata=zeros(num_CL,6);
new_crosslinkdata(:,2)=2; new_crosslinkdata(:,3)=2;
new_crosslinkdata(:,4)= xlow +(xup-xlow).*rand(num_CL,1);  new_crosslinkdata(:,5)= ylow +(yup-ylow).*rand(num_CL,1);   new_crosslinkdata(:,6)= zlow +(zup-zlow).*rand(num_CL,1);

cal1 = sqrt(new_crosslinkdata(:,4).^2+new_crosslinkdata(:,5).^2 + new_crosslinkdata(:,6).^2  );
todelete = find(cal1>(clot_diameter/2));          new_crosslinkdata(todelete,:)=[];
new_crosslinkdata(:,1)=1:length(new_crosslinkdata(:,1));
crosslinkData=new_crosslinkdata;
num_CL=length(crosslinkData(:,1));

% calculating # of plt , and generating location
num_PLT= ceil(num_p_ref * cubicV /(clotsize_ref^3));  
xup = clot_sidelength /2; xlow =- clot_sidelength/2;   yup = clot_sidelength /2; ylow = -clot_sidelength/2;    zup = clot_sidelength/2;  zlow = -clot_sidelength/2;

pdata=zeros(num_PLT,5);
pdata(:,1)=3;   pdata(:,2)=3;
pdata(:,3)= xlow +(xup-xlow).*rand(num_PLT,1);   pdata(:,4)= ylow +(yup-ylow).*rand(num_PLT,1);      pdata(:,5)= zlow +(zup-zlow).*rand(num_PLT,1);
cal2 = sqrt(pdata(:,3).^2+pdata(:,4).^2+ pdata(:,5).^2);
todelete = find(cal2>(clot_diameter/2));   pdata(todelete,:)=[];
num_PLT=length(pdata(:,1));

pdata(:,3)=pdata(:,3).*0.9; pdata(:,4)=pdata(:,4).*0.9; 
% add in PLT body beads
load('bonds61');        
load('data_coords61');  
data_coords_p=data_coords61;     bonds_p=bonds61;
ss=0.2;
data_coords_p(:,4:6)=data_coords_p(:,4:6).*ss;
BLcheck1 = (data_coords_p(1,4)-data_coords_p(5,4))^2+(data_coords_p(1,5)-data_coords_p(5,5))^2+(data_coords_p(1,6)-data_coords_p(5,6))^2;
BLcheck2 = (data_coords_p(58,4)-data_coords_p(59,4))^2+(data_coords_p(58,5)-data_coords_p(59,5))^2+(data_coords_p(58,6)-data_coords_p(59,6))^2;
plt_BL1=sqrt(BLcheck1);  plt_BL2=sqrt(BLcheck2);

plateletdata = pdata;   pdata = [];
center_id = []; center_id = find(data_coords_p(:,4)==0 & data_coords_p(:,5)==0 & data_coords_p(:,6)==0);
for i = 1:num_PLT
    add_pdata = zeros(length(data_coords_p),5);
    add_pdata(:,1)=3;  add_pdata(:,2)=3;
    add_pdata(:,3)= plateletdata(i,3)+data_coords_p(:,4);
    add_pdata(:,4)= plateletdata(i,4)+data_coords_p(:,5);
    add_pdata(:,5)= plateletdata(i,5)+data_coords_p(:,6);
    add_pdata(center_id,1)=0;  add_pdata(center_id,2)=0;
    pdata = [pdata;add_pdata];
end

pc = find(pdata(:,1)==0);
type_vec = zeros(num_PLT,1);
k = ceil(num_PLT/20);
rest = num_PLT-(k-1)*20;
for i=1:k-1
    iii =( 1+(i-1)*20 ) : (20+(i-1)*20 );
    type_vec(iii) = 4:23;
end
iii = (1+(k-1)*20 ):num_PLT;
type_vec(iii) = 4:(4+length(iii)-1);
pdata(pc,1)=type_vec;   pdata(pc,2)=type_vec;
sprintf('Num PLT in clot=%d, plt bond L=%.3f&%.3f',num_PLT,plt_BL1,plt_BL2)



kk = find(pdata(:,2)>=4 &pdata(:,2)<=53 );


%%
[ NatomN, MatomN, TmolN, TatomN, CxN, CyN, CzN, TbondN, N1bondN, N2bondN, TangN, N1angN, N2angN, N3angN, Tdih,N1dih,N2dih,N3dih,N4dih, Spacing_nodes, FilamentDistribution] ...
    = membrane_gel_folding_sphere( Smap, Box, Initial_Network_Rho, Leq, ave_connectivity, r_cutoff, H, r1, MatomOne, periodic_x, periodic_y, periodic_z, ...
    atom_type, mol_type, bond_type, angle_type, psize, gelsize, percentActive, ...
    pdata, crosslinkData,p_cutoff,maxConnections,num_PLT,r_cutoff_inner,clot_sidelength,num_rbc,box_sidelength,data_coords_p,bonds_p);

[Natom Matom Tmol Tatom Cx Cy, Cz Tbond N1bond N2bond Tang N1ang N2ang N3ang]...
    = add_arrays( Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang,...
    NatomN, MatomN, TmolN, TatomN, CxN, CyN, CzN, TbondN, N1bondN, N2bondN, TangN, N1angN, N2angN, N3angN);

%%

FileNameNew = sprintf('clot_%d_%d_%g_%g_%g.dat',Leq, ave_connectivity, Initial_Network_Rho,r_cutoff,p_cutoff);

create_output_file_new(FileNameNew, Box, Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang,Tdih,N1dih,N2dih,N3dih,N4dih,box_sidelength);

