function [] = Graph_UB_sol(xx, P, polyL, Int_i, teta, Origine, scale, varargin)
%GRAPH_UB_SOL Summary of this function goes here
%   Detailed explanation goes here
%%
if isempty(varargin)
    figure(2),box on, axis equal, hold on % put the figure in a function
    ax=gca();
    set(gca,'Box','on');
    xticklabels({});
    yticklabels({});
else
    ax=varargin{1};
end
% title(ax,['Deformed shape (\lambda=', num2str(xx.pval,'%.4f'), ' ,\theta=', num2str(teta*180/pi,'%.0f'), '°)'])
% xlabel(ax,'x [mm]'),ylabel(ax,'y [mm]')

hold(ax,"on")
% axis(ax,"equal")

Nvertmax=size(polyL,2)-6; % delete columns with xg yg zg Area thickness density
% for i=1:size(polyL,1)
%     for k=1:Nvertmax
%         if polyL(i,k)~=0
%             coord_x(k)=P(polyL(i,k),2);    
%             coord_y(k)=P(polyL(i,k),3);            
%         end           
%     end 
%     polyline=polyshape(coord_x,coord_y);
%     plot(polyline,'FaceColor','black','FaceAlpha',0.1),hold on
%     clear coord_x
%     clear coord_y
% end

% Plot deformed configuration
U_phi(1:size(polyL,1)*3,1)=xx.ub(1:size(polyL,1)*3,1);
U_phi=U_phi*scale;
colorbar(ax);
for r=1:size(polyL,1)
    Uvec=U_phi(3*(r-1)+1:3*r,1);
    Ux=Uvec(1);
    Uy=Uvec(2);
    phi=Uvec(3);
    G=polyL(r,Nvertmax+1:Nvertmax+3);
    cf=[];
    for i=1:Nvertmax
        Nk=polyL(r,i);
        if polyL(r,i)~=0
           Pk=P(Nk,2:4); 
           DU=[Ux Uy 0]+([0 -phi 0;phi 0 0;0 0 0]*(Pk-G)')';
           Pnew(i,1:2)=P(polyL(r,i),2:3)+1*DU(1:2); %Pnew deformed node
           R=[cos(teta), sin(teta);
             -sin(teta), cos(teta)];
           P_trasl=Pnew(i,1:2) - Origine;    
           P_rot=(R * P_trasl')';            
           Pnew(i,1:2)=P_rot+Origine; %Pnew deformed and rotated node
           cf=[cf;norm(DU(1:2))];
        end 
    end

    % Pnew=[Pnew; Pnew(1,:)];
    % fill(ax,Pnew(:,1),Pnew(:,2), [240/255 240/255 240/255])
    patch(ax,Pnew(:,1),Pnew(:,2),cf,'FaceColor','interp')
    clear Pnew
end
% - add for adaption of old MATLAB version
pb=get(ax,"PlotBoxAspectRatio");
set(ax,"PlotBoxAspectRatioMode","manual");
set(ax,"PlotBoxAspectRatio",pb);
set(ax,"DataAspectRatio",[1 1 1]);
end

function Aeq=eval_Aeq(P,polyL,Int_i)

Gass=zeros(4*size(Int_i,1),3*size(polyL,1));
Mass=zeros(4*size(Int_i,1),4*size(Int_i,1));

for i=1:size(Int_i,1)
    % number of nodes and number of elements
    N1=P(Int_i(i,1));
    N2=P(Int_i(i,2));
    E1=Int_i(i,3);
    E2=Int_i(i,4);
    ni=Int_i(i,8:10);
    ti=Int_i(i,5:7);

    Gk_n1_E1=[1 0 -(P(N1,3)-polyL(E1,size(polyL,2)-4))
              0 1 (P(N1,2)-polyL(E1,size(polyL,2)-5))];

    Gk_n2_E1=[1 0 -(P(N2,3)-polyL(E1,size(polyL,2)-4))
              0 1 (P(N2,2)-polyL(E1,size(polyL,2)-5))];

    Gk_n1_E2=[1 0 -(P(N1,3)-polyL(E2,size(polyL,2)-4))
              0 1 (P(N1,2)-polyL(E2,size(polyL,2)-5))];

    Gk_n2_E2=[1 0 -(P(N2,3)-polyL(E2,size(polyL,2)-4))
              0 1 (P(N2,2)-polyL(E2,size(polyL,2)-5))];

    Gloc_E1=[ni(1:2)*Gk_n1_E1
             ti(1:2)*Gk_n1_E1
             ni(1:2)*Gk_n2_E1
             ti(1:2)*Gk_n2_E1];

    Gloc_E2=[-ni(1:2)*Gk_n1_E2
             -ti(1:2)*Gk_n1_E2
             -ni(1:2)*Gk_n2_E2
             -ti(1:2)*Gk_n2_E2];

    Mk=[-Int_i(i,13) -Int_i(i,13) 0 0
        -1 1 0 0
        0 0 -Int_i(i,13) -Int_i(i,13)
        0 0 -1 1];

    k=i;    
    j=E1;
    Gass(4*(k-1)+1:4*k,3*(j-1)+1:3*j)=Gloc_E1;
    j=E2;
    Gass(4*(k-1)+1:4*k,3*(j-1)+1:3*j)=Gloc_E2;

    Mass(4*(k-1)+1:4*k,4*(k-1)+1:4*k)=Mk;

end

Aeq=[Gass Mass];
end

