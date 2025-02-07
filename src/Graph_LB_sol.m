function [] = Graph_LB_sol(xx, P, polyL, Int_i, teta, Origine, scale, varargin)
%GRAPH_UB_SOL Summary of this function goes here
%   Detailed explanation goes here
%%
if isempty(varargin)
    figure('Name','Force'),box on, axis equal, hold on %mettere la figura in una funzione
    set(gca,'Box','on');
    xticklabels({});
    yticklabels({});
    daspect([1 1 1]);
    ax=gca();
else
    ax=varargin{1};
end

hold(ax,"on");
% axis(ax,"equal");


% plot resultant
xf=reshape(xx.lb(2:(end-3)),4,[])';
R=[cos(teta), sin(teta);
  -sin(teta), cos(teta)];
charlength=sqrt(mean(polyL(:,end-2)));
for i=1:size(polyL,1)
    poly=polyL(i,1:end-6);
    poly=poly(poly~=0);
    nds=P(poly,2:3);             
    nds=(R * (nds - Origine)')'+Origine;
    fill(ax,[nds(:,1);nds(1,1)]',[nds(:,2);nds(1,2)]', [240/255 240/255 240/255])
    % line(ax,[nds(:,1);nds(1,1)]',[nds(:,2);nds(1,2)]','color','k','linewidth',0.5);
end
colorbar(ax);
for i=1:size(Int_i,1)
    Int_ctd=(P(Int_i(i,1),2:3)+P(Int_i(i,2),2:3))/2;
    Int_len=Int_i(i,end-2);
    Int_vecs=Int_i(i,5:6);
    Int_vecn=Int_i(i,8:9);
    fn=xf(i,1)+xf(i,3);
    fs=xf(i,2)+xf(i,4);
    cf=norm(fn*Int_vecn+fs*Int_vecs,2);
    if cf>max(max(abs(xf)))*1e-9
        % ff=log10(cf)*normalize(fn*Int_vecn+fs*Int_vecs,2,"norm");
        ff=fn*Int_vecn+fs*Int_vecs;
        ee=(-xf(i,1)*Int_len/2+xf(i,3)*Int_len/2)/fn;
        xs=Int_ctd(1)+ee*Int_vecs(1)-scale*charlength*ff(1)/2;
        ys=Int_ctd(2)+ee*Int_vecs(2)-scale*charlength*ff(2)/2;
        xe=Int_ctd(1)+ee*Int_vecs(1)+scale*charlength*ff(1)/2;
        ye=Int_ctd(2)+ee*Int_vecs(2)+scale*charlength*ff(2)/2;
        % pts=[xs ys];
        % pte=[xe ye];
        pts=(R * ([xs ys] - Origine)')'+Origine;
        pte=(R * ([xe ye] - Origine)')'+Origine;
        % patch(ax,[pts(1,1),pte(1,1)],[pts(1,2),pte(1,2)],cf,'EdgeColor','flat','LineWidth',1.5);
        surface(ax,[pts(1,1) pte(1,1);pts(1,1) pte(1,1)],[pts(1,2) pte(1,2); pts(1,2) pte(1,2)],[0 0;0 0],[cf cf;cf cf],'facecol','none','edgecol','interp','linew',1.5);
    end
end
% title(ax,['Force line (\lambda=', num2str(-xx.dval,'%.4f'), ' ,\theta=', num2str(teta*180/pi,'%.0f'), '°)'])
% title(sprintf('\\lambda=%.4f',-fval));
% - add for adaption of old MATLAB version
pb=get(ax,"PlotBoxAspectRatio");
set(ax,"PlotBoxAspectRatioMode","manual");
set(ax,"PlotBoxAspectRatio",pb);
set(ax,"DataAspectRatio",[1 1 1]);
end

