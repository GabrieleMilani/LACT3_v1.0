%-------------------------------------------------------------------------%
%           LIMIT ANALYSIS CODE FOR TILTING TABLE TESTS (LACT3)           %
%-------------------------------------------------------------------------%
% This MATLAB script package, featuring a user-interactive GUI, produces  % 
% a tool for the fast computation of tilting table tests of masonry       %
% systems based on limit analysis.  Check the "README.md" file for more   %
% instructions and details.                                               %
%-------------------------------------------------------------------------%
% Developed by:                                                           %
%                          Yiwei Hua,                                     %
%                          Martina Buzzetti,                              % 
%                          Natalia Pingaro,                               %
%                          Luis C.M. da Silva,                            %
%                          Gabriele Milani                                %
%                                                           April 2, 2025 %
%-------------------------------------------------------------------------%
% Cited as:                                                               %
% Y. Hua, M. Buzzetti, N. Pingaro, L.C.M. da Silva, G. Milani. A compu-   %
% terized tool for the kinematic limit analysis of 2D masonry structures  %
% failing on a tilting table.                                             %
%-------------------------------------------------------------------------%

%% ------------------------- INITIALIZING GUI -------------------------- %%
% function GUI_EXAMPLE
clear all
clc;
warning off
all_fig = findall(0, 'type', 'figure');
close(all_fig)
file = mfilename('fullpath');
[current_path,name]=fileparts(file);
cd(current_path);
addpath(genpath(current_path));
%% Configuration setting
% figure
global lb1 ax1 ax2 ax3 ax4 ef1 ef2 ef3 ef4 ef5 ef6 dd tg t4 tt1 tt2 sp fig
ss=get(0,"ScreenSize");
fig = uifigure('Position',[150,150,ss(3:4)*0.75],"Name","Tilt test analysis");
gl = uigridlayout(fig,[3 3]);
gl.ColumnWidth={".8x","1x","1x"};
gl.RowHeight={"1.2x","0.8x","0.1x"};
% panel
p1 = uipanel(gl);
p1.Layout.Row=1;
p1.Layout.Column=1;
p1.Title="Model plot";
g1 = uigridlayout(p1,[1 1]);
ax1=uiaxes(g1,"XTickLabel",{},"YTickLabel",{},"XTick",{},"YTick",{},"Box","on");
p2 = uipanel(gl);
p2.Layout.Row=1;
p2.Layout.Column=2;
p2.Title="Upper bound solution: Velocity field";
g2 = uigridlayout(p2,[1 1]);
ax2=uiaxes(g2,"XTickLabel",{},"YTickLabel",{},"XTick",{},"YTick",{},"Box","on");
p3 = uipanel(gl);
p3.Layout.Row=1;
p3.Layout.Column=3;
p3.Title="Lower bound solution: Force line plot";
g3 = uigridlayout(p3,[1 1]);
ax3=uiaxes(g3,"XTickLabel",{},"YTickLabel",{},"XTick",{},"YTick",{},"Box","on");
ax4 = uiaxes(gl,"Box","on");
ax4.XLabel.String="Tilt angle \it\theta";
ax4.XLabel.FontSize=13;
ax4.YLabel.String="Collaspe multiplier \alpha";
ax4.YLabel.FontSize=13;
ax4.Layout.Row=[2 3];
ax4.Layout.Column=[2 3];
grid(ax4,"on");

tg = uitabgroup(gl);
tg.Layout.Row=2;
tg.Layout.Column=1;
t1 = uitab(tg,"Title","Model input");
t2 = uitab(tg,"Title","Material settings");
t3 = uitab(tg,"Title","Iteration settings");
t4 = uitab(tg,"Title","Output messages");

b1 = uibutton (gl,"Text","Run");
b1.Layout.Row=3;
b1.Layout.Column=1;
b1.ButtonPushedFcn=@RunButtonPushed;

g4 = uigridlayout(t1,[3 2]);
g4.RowHeight={"fit","fit","1x"};
b2 = uibutton (g4,"Text","Select file");
b2.Layout.Row=1;
b2.Layout.Column=1;
b2.ButtonPushedFcn=@SleButtonPushed;
b3 = uibutton (g4,"Text","Plot model");
b3.Layout.Row=1;
b3.Layout.Column=2;
b3.ButtonPushedFcn=@PltButtonPushed;
% lb1=uilabel(g4,"Text",sprintf(">> File path:\n"),"WordWrap","on");
lb1=uilabel(g4,"Text",sprintf(">> File path:\n"));
lb1.Layout.Row=2;
lb1.Layout.Column=[1 2];
tt2=uitextarea(g4,"Enable","off");
tt2.Layout.Row=3;
tt2.Layout.Column=[1 2];

g5 = uigridlayout(t2,[4 2]);
g5.RowHeight={"fit","fit","fit","fit"};
lb2=uilabel(g5,"Text","Friction angle");
lb2.Layout.Row=1;
lb2.Layout.Column=1;
lb3=uilabel(g5,"Text","Cohesion");
lb3.Layout.Row=2;
lb3.Layout.Column=1;
lb9=uilabel(g5,"Text","Density");
lb9.Layout.Row=3;
lb9.Layout.Column=1;
ef1 = uieditfield(g5,"numeric","ValueDisplayFormat","%.2f\x00B0","Limits",[0 Inf],"Value",30);
ef1.Layout.Row=1;
ef1.Layout.Column=2;
ef2 = uieditfield(g5,'numeric',"ValueDisplayFormat","%.2fMPa","Limits",[0 Inf]);
ef2.Layout.Row=2;
ef2.Layout.Column=2;
ef6 = uieditfield(g5,"numeric","ValueDisplayFormat","%.2fkg/m\x00B3","Limits",[0 Inf],"Value",2500);
% set(ef6, 'Interpreter', 'tex');
ef6.Layout.Row=3;
ef6.Layout.Column=2;

g6 = uigridlayout(t3,[4 2]);
g6.RowHeight={"fit","fit","fit","fit"};
lb4=uilabel(g6,"Text","Initial tilt angle");
lb4.Layout.Row=1;
lb4.Layout.Column=1;
lb5=uilabel(g6,"Text","Initial step increment");
lb5.Layout.Row=2;
lb5.Layout.Column=1;
lb6=uilabel(g6,"Text","Converging tolerance");
lb6.Layout.Row=3;
lb6.Layout.Column=1;
lb7=uilabel(g6,"Text","Step control");
lb7.Layout.Row=4;
lb7.Layout.Column=1;
ef3 = uieditfield(g6,"numeric","ValueDisplayFormat","%.2f\x00B0","Limits",[0 Inf]);
ef3.Layout.Row=1;
ef3.Layout.Column=2;
ef4 = uieditfield(g6,'numeric',"ValueDisplayFormat","%.2f\x00B0","Limits",[1e-3 Inf],"Value",1);
ef4.Layout.Row=2;
ef4.Layout.Column=2;
ef5 = uieditfield(g6,"numeric","ValueDisplayFormat","%.2e","Limits",[1e-6 Inf],"Value",0.01);
ef5.Layout.Row=3;
ef5.Layout.Column=2;
dd = uidropdown(g6,"Items",["Automatic","Fixed"]);
dd.Layout.Row=4;
dd.Layout.Column=2;

g7 = uigridlayout(t4,[2 2]);
g7.RowHeight={"1x","fit"};
tt1=uitextarea(g7,"Enable","off");
tt1.Layout.Row=1;
tt1.Layout.Column=[1 2];
sp=uispinner(g7,"Enable","off","Limits",[1 10]);
sp.Layout.Row=2;
sp.Layout.Column=2;
sp.ValueChangingFcn=@SpiValChanged;
lb8=uilabel(g7,"Text","Iteration frame");
lb8.Layout.Row=2;
lb8.Layout.Column=1;

function SleButtonPushed(src,event)
    global lb1 filename pathname
    [filename, pathname]=SelectFile(lb1);
end

function PltButtonPushed(src,event)
    global filename pathname P polyL int_i ax1 tt1 tt2 ax2 ax3 ax4 sp
    cla(ax1);
    % axis(ax1,"normal")
    axis(ax1,"auto")
    [P, polyL, int_i] = Load_wall_data_2(filename, pathname, ax1);%
    % axis(ax1,"equal")
    tt2.Enable="on";
    % Write output message
    tt2.Value=sprintf(join(["------------------Model geometry summary------------------",...
                            ">> Software version: v1.0",...
                            ">> Problem dimension: 2D",...
                            ">> Number of nodes: %d",...
                            ">> Number of elements: %d",...
                            ">> Number of interfaces: %d"],'\n'),...
                            size(P,1), size(polyL,1),size(int_i,1));
    % Graphic plot
    cla(ax2);
    axis(ax2,'auto')
    cla(ax3);
    axis(ax3,'auto')
    cla(ax4);
    tt1.Value="";
    tt1.Enable="off";
    sp.Value=1;
    sp.Enable="off";
end

function RunButtonPushed(src,event) 
    global ef1 ef2 ef3 ef4 ef5 ef6 dd ax2 ax3 ax4 P polyL int_i tg t4 tt1 sp xx Converged_step fig scale_u scale_f
    xx=struct();
    d=uiprogressdlg(fig,'Title','Run analysis','Indeterminate','on','Cancelable','on');
    d.Message = 'Initial setting';
    drawnow
    % get para
    phi = ef1.Value;
    c = ef2.Value;
    tol = ef5.Value*pi/180;
    step_ctl = dd.Value;
    cla(ax2);
    axis(ax2,'auto')
    cla(ax3);
    axis(ax3,'auto')
    cla(ax4);
    % axis(ax2,"normal")
    % axis(ax3,"normal")
    % Material assignment
    [polyL, int_i] = Material(polyL, int_i, c, phi, 1, ef6.Value);
    MaxIterNum=90;
    step = ef4.Value*pi/180;
    xx(1).teta = ef3.Value*pi/180;
    i=1;
    Last_converged_teta=0;
    while i<MaxIterNum && xx(i).teta<pi/2
        d.Message = sprintf('Solving iteration %d',i);
        [xx(i).ub,xx(i).pval] = Solver_in_UB(P, polyL, int_i, xx(i).teta);
        [xx(i).lb,xx(i).dval] = Solver_in_LB(P, polyL, int_i, xx(i).teta);
        if xx(i).pval>0
            xx(i).is_Converged=true;
            Last_converged_teta=xx(i).teta;
            % d.Message = sprintf('Solving iteration %d...converged',i);
        else
            xx(i).is_Converged=false;
            % d.Message = sprintf('Solving iteration %d...unconverged',i);
        end
        % compute step size
        if step>tol
            i=i+1;
            step = Step_control(step,xx,step_ctl);
            xx(i).teta = min(Last_converged_teta + step, pi/2);            
        else
            break;
        end
        if d.CancelRequested
            break;
        end
    end
    Converged_step=find([xx.is_Converged]);
    if ~isempty(Converged_step)
        d.Message = 'Converged solution found';
        pause(.5);
        d.Message = 'Plot results';
        frame=Converged_step(end);
        % components update
        tg.SelectedTab=t4;
        sp.Enable="on";
        sp.Limits=[1 length(Converged_step)];
        sp.Value=length(Converged_step);
        sp.ValueChangingFcn=@SpiValChanged;
        tt1.Enable="on";
        % Write output message
        tt1.Value=sprintf(join(["------------------Collapse results summary------------------",...
                                ">> Limit solution:",...
                                "   - Collapse tilt angle: %.2f\x00B0",...
                                "   - Upper bound load: %.4f",...
                                "   - Lower bound load: %.4f",...
                                ">> Current frame: %d",...
                                "   - tilt angle: %.2f\x00B0",...
                                "   - Upper bound load: %.4f",...
                                "   - Lower bound load: %.4f"],'\n'),...
                                xx(frame).teta/pi*180, xx(frame).pval,-xx(frame).dval,...
                                frame,xx(frame).teta/pi*180, xx(frame).pval,-xx(frame).dval);
        % Graphic plot
        aabb_x=max(P(:,2))-min(P(:,2));
        aabb_y=max(P(:,3))-max(P(:,3));
        scale_u=0.1*max(aabb_x,aabb_y)/max(abs(xx(frame).ub));
        scale_f=0.5*max(aabb_x,aabb_y)/max(abs(xx(frame).lb));
        % [~,Origine]=eval_BC(P,polyL,int_i);
        Origine=Get_Rotation_origin(P);
        Graph_UB_sol(xx(frame), P, polyL, int_i, xx(frame).teta, Origine, scale_u, ax2)
        Graph_LB_sol(xx(frame), P, polyL, int_i, xx(frame).teta, Origine, scale_f, ax3)
        % axis(ax2,"equal")
        % axis(ax3,"equal")
        % Iteration plot
        plot(ax4,[xx(Converged_step).teta]'/pi*180,[xx(Converged_step).pval]','.-',"LineWidth",1,"MarkerSize",10);
        xline(ax4,xx(frame).teta/pi*180,'-',{'Current frame'},"Color",[0.6350 0.0780 0.1840],"LineWidth",1.5)
        % xlim(ax4,"tight");
        set(ax4,"XLim",[xx(Converged_step(1)).teta/pi*180 xx(Converged_step(end)).teta/pi*180+1]);
        close(d);
    else
        close(d);
        uialert(fig,"Initial configuration is not feasible in static. Please change the setting","No feasible solution");
    end
end

function SpiValChanged(src,event)
    global ax2 ax3 ax4 P polyL int_i xx Converged_step tt1 scale_u scale_f
    frame=Converged_step(event.Value);
    % Write output message
    tt1.Value=sprintf(join(["------------------Collapse results summary------------------",...
                            ">> Limit solution:",...
                            "   - Collapse tilt angle: %.2f\x00B0",...
                            "   - Upper bound load: %.4f",...
                            "   - Lower bound load: %.4f",...
                            ">> Current frame: %d",...
                            "   - tilt angle: %.2f\x00B0",...
                            "   - Upper bound load: %.4f",...
                            "   - Lower bound load: %.4f"],'\n'),...
                            xx(Converged_step(end)).teta/pi*180, xx(Converged_step(end)).pval,-xx(Converged_step(end)).dval,...
                            frame,xx(frame).teta/pi*180, xx(frame).pval,-xx(frame).dval);
    % Graphic plot 
    % [~,Origine]=eval_BC(P,polyL,int_i);
    Origine=Get_Rotation_origin(P);
    cla(ax2);
    axis(ax2,'auto')
    cla(ax3);
    axis(ax3,'auto')
    Graph_UB_sol(xx(frame), P, polyL, int_i, xx(frame).teta, Origine, scale_u, ax2)
    Graph_LB_sol(xx(frame), P, polyL, int_i, xx(frame).teta, Origine, scale_f, ax3)
    % axis(ax2,"equal")
    % axis(ax3,"equal")
    % Iteration plot
    delete(ax4.Children(1));
    xline(ax4,xx(frame).teta/pi*180,'-',{'Current frame'},"Color",[0.6350 0.0780 0.1840],"LineWidth",1.5);
end

function [Origine]=Get_Rotation_origin(P)
    P_base=P(P(:,3)>=-1e-3 & P(:,3)<=1e-3,2:3);
    % Find and save tilting plane origin
    [~, maxIndex] = max(P_base(:,1));
    Origine = P_base(maxIndex, :);
end