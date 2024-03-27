
function [bondListPlatelets, final_ghostdata] =LinkPlateletsAndNetwork(output1111)

Leeeq = 0.045;
platelets = output1111(output1111(:,2) == 3, :);
crosslinks = output1111(output1111(:,2) == 2, :);
maxConnections = 5;
cutoff = 4.5;

bondListPlatelets = [];
local_connections = [];

ghost_data=[];  %[ghostID, x,y,z]
local_ghost=[];
ghostID= length(output1111);

local_connections_PC = [];
local_connections_PG = [];
local_connections_GG = [];
local_connections_GC = [];


bondList_PC= [];
bondList_PG= [];
bondList_GG= [];
bondList_GC= [];

ghostforthisbond = [];



for i = 1:2
    dist = sqrt((platelets(i,3) - crosslinks(:,3)).^2 + (platelets(i,4) - crosslinks(:,4)).^2 + (platelets(i,5) - crosslinks(:,5)).^2);
    ind = find(dist < cutoff);
%     m=length(ind);
%     choose = randperm(m);
%     for i = 1:m
%         ind(i) = ind(choose(i));
%     end
    
    for j = 1:length(ind)  % loop for all found close crosslinks
        if length(local_connections) < maxConnections
            local_connections = [local_connections; platelets(i,1), crosslinks(ind(j))]; 
            d=sqrt((platelets(i,3) - crosslinks(ind(j),3)).^2 + (platelets(i,4) - crosslinks(ind(j),4)).^2 + (platelets(i,5) - crosslinks(ind(j),5)).^2);
            if d>=Leeeq   % need to seperate the bonds 
                NumGhost = round( d / Leeeq )-1;
                for k=1:NumGhost
                    ghostID=ghostID+1;
                    gx = platelets(i,3)+ (k/(NumGhost+1))*(crosslinks(ind(j),3)-platelets(i,3));
                    gy = platelets(i,4)+ (k/(NumGhost+1))*(crosslinks(ind(j),4)-platelets(i,4));
                    gz = platelets(i,3)+ (k/(NumGhost+1))*(crosslinks(ind(j),5)-platelets(i,5));
                    ghostforthisbond = [ghostforthisbond; ghostID, gx,gy,gz ];
                    local_ghost = [local_ghost; ghostID,gx,gy,gz];
                end  
                local_connections_PG = [local_connections_PG; platelets(i,1), ghostforthisbond(1,1)];
                local_connections_GC = [local_connections_GC; ghostforthisbond(NumGhost,1), crosslinks(ind(j),1)];
                for m=1:NumGhost-1
                    local_connections_GG =[local_connections_GG; ghostforthisbond(m,1),ghostforthisbond(m+1,1)];
                end  
            elseif d<Leeeq
                local_connections_PC = [local_connections_PC; platelets(i,1), crosslinks(ind(j),1)];
            end  
        end
        ghostforthisbond = [];

    end
    bondListPlatelets = [bondListPlatelets; local_connections];
    
    bondList_PC= [bondList_PC; local_connections_PC];
    bondList_PG= [bondList_PG; local_connections_PG];
    bondList_GG= [bondList_GG; local_connections_GG];
    bondList_GC= [bondList_GC; local_connections_GC];
    ghost_data = [ghost_data; local_ghost];

    
    local_connections = [];
    local_ghost=[];
    local_connections_PC=[];
    local_connections_PG = [];
    local_connections_GG = [];
    local_connections_GC = [];

end


bondListPlatelets = [bondList_PC;bondList_PG;bondList_GG;bondList_GC];
bondListPlatelets= sort(bondListPlatelets,2 );
bondListPlatelets = unique(bondListPlatelets,'rows');


   %   ghost_data=[];  %[ghostID, x,y,z]  
   ng=length(ghost_data);
   type = 3*ones(ng,1);
   final_ghostdata=zeros(ng,6);
   final_ghostdata(:,1) = ghost_data(:,1);
   final_ghostdata(:,2) =type;
   final_ghostdata(:,3) =ghost_data(:,2);
   final_ghostdata(:,4) =ghost_data(:,3);
   final_ghostdata(:,5) =ghost_data(:,4);
   final_ghostdata(:,6) = type;
