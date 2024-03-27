function[my2,my1]=prepare20220113(output1111,output22222,H,gelsize,p_cutoff,maxConnections,num_per_p,num_p)

L.x = gelsize;
L.y = H;
L.z = gelsize;
ccc = find(output1111(:,2) >= 4 & output1111(:,2) <=23);
p_centers = output1111(ccc,: );

if length(p_centers(:,1))~= num_p
    disp('error: wrong plt center finder')
end
crosslinks = output1111(output1111(:,2) == 2, :);

bondtype = 3;

numberofPlatelet = num_p;
bondID = length(output22222);
localconnection = [];
bondcount = 0;
my2=[];

plateletmol_starter = max(output1111(:,6));
mypos = find(output1111(:,2)>=3  );

numberperplatelet = num_per_p;

add_bonds_tokeep=[];

for n = 1:numberofPlatelet
    plateletmol_starter = plateletmol_starter+1;
    
    %% change mol of pl
    changeid = mypos(1:numberperplatelet);
    output1111(changeid,6) = plateletmol_starter;
    mypos(1:numberperplatelet)=[];
    
    %% start to give bonds between platelet and fibrin
    usedfilament =[];
    thisplatelet = p_centers(n,:);
    
    for j = 1:length(crosslinks(:,1))
        dist(j) = sqrt ( ( min(abs(thisplatelet(1,3) - crosslinks(j,3)),L.x-abs(thisplatelet(1,3) - crosslinks(j,3))))^2 ...
            + ( min(abs(thisplatelet(1,4) - crosslinks(j,4)),L.y-abs(thisplatelet(1,4) - crosslinks(j,4))))^2 ...
            + ( min(abs(thisplatelet(1,5) - crosslinks(j,5)),L.z-abs(thisplatelet(1,5) - crosslinks(j,5))))^2);
    end
    
    % add a bond to keep the platelet not escaped
     [mindist(n),where]=min(dist);
     bondID = bondID+1;
    localconnection = [localconnection; bondID,3,thisplatelet(1,1),crosslinks(where,1),mindist(n)];
 
 
    
end
my2 = [my2;localconnection];
my1 = output1111;


end