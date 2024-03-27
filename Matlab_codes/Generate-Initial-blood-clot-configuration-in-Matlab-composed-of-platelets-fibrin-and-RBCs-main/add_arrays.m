function [ Natom Matom Tmol Tatom Cx Cy Cz Tbond N1bond N2bond Tang N1ang N2ang N3ang]...
    = add_arrays( Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang,...
            NatomS, MatomS, TmolS, TatomS, CxS, CyS, CzS, TbondS, N1bondS, N2bondS, TangS, N1angS, N2angS, N3angS)

% [msx msy msz] = size(Map);
% if(msx == 1 || msy == 1 || msz == 1)
%     disp(sprintf('Map[%d, %d, %d] size error', msx, msy, msz));
% end

Natom_cur = length(Natom);
Natom(Natom_cur+1:Natom_cur+length(NatomS),1) = Natom_cur + NatomS;
Matom(Natom_cur+1:Natom_cur+length(NatomS),1) = MatomS;
Tmol(Natom_cur+1:Natom_cur+length(NatomS),1) = TmolS;
Tatom(Natom_cur+1:Natom_cur+length(NatomS),1) = TatomS;

Cx(length(Cx)+1:length(Cx)+length(CxS),1) = CxS;
Cy(length(Cy)+1:length(Cy)+length(CyS),1) = CyS;
Cz(length(Cz)+1:length(Cz)+length(CzS),1) = CzS;

Nbond_cur = length(Tbond);
Tbond(Nbond_cur+1:Nbond_cur+length(TbondS),1) = TbondS;
N1bond(Nbond_cur+1:Nbond_cur+length(TbondS),1) = Natom_cur + N1bondS;
N2bond(Nbond_cur+1:Nbond_cur+length(TbondS),1) = Natom_cur + N2bondS;

Nang_cur = length(Tang);
Tang(Nang_cur+1:Nang_cur+length(TangS),1) = TangS;
N1ang(Nang_cur+1:Nang_cur+length(TangS),1) = Natom_cur + N1angS;
N2ang(Nang_cur+1:Nang_cur+length(TangS),1) = Natom_cur + N2angS;
N3ang(Nang_cur+1:Nang_cur+length(TangS),1) = Natom_cur + N3angS;



