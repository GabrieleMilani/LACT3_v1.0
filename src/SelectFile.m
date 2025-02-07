function [filename, pathname]=SelectFile(lb1)
% global lb1 filename pathname

    [filename, pathname] = uigetfile({'*.dxf'},'Open File','Multiselect','off'); % choose file to open
        if filename==0
            errordlg('no file selected','file error');return
        end
    % lb1.Text=sprintf('File path:\n%s%s',pathname,filename);
    % pathname=[pathname(1:30),'\n',pathname(31:end)];
    lb1.Text=sprintf('>> File path: %s\n     %s%s',pathname(1:40),pathname(41:end),filename);
end