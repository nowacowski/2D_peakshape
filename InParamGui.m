function new_params = InParamGui(params)
    param_fit = reshape(params,6,[])';
    cancel_change = param_fit;

% Main window
    hwin = figure('Visible','off');
    hwin.NumberTitle = 'off';
    hwin.Name = 'Input initial parameters';
    hwin.MenuBar = 'none';
    hwin.Resize = 'off';
    hwin.Units = 'normalized';
    hwin.Position = [0.1 0.1 0.4 0.3];

    bgclr = get(hwin,'Color');

% Table
    htable = uitable(hwin);
    htable.Data = param_fit;
    htable.ColumnName = {'A', 'mu_x', 'w_x', 'mu_y', 'w_y', 'rho'};
    htable.ColumnEditable = true;
    htable.Units = 'normalized';
    htable.Position = [0.05 0.15 0.9 0.8];

% Add row
    hrow = uicontrol('Style', 'pushbutton');
    hrow.BackgroundColor = bgclr;
    hrow.String = 'Add row';
    hrow.Callback = @RowFcn;
    hrow.Units = 'normalized';
    hrow.Position = [0.2 0.03 0.15 0.1];
    
% remove row
    hrowrem = uicontrol('Style', 'pushbutton');
    hrowrem.BackgroundColor = bgclr;
    hrowrem.String = 'Remove row';
    hrowrem.Callback = @RowRemFcn;
    hrowrem.Units = 'normalized';
    hrowrem.Position = [0.4 0.03 0.15 0.1];

% Cancel
    hcancel = uicontrol('Style','pushbutton');
    hcancel.BackgroundColor = bgclr;
    hcancel.String = 'Cancel';
    hcancel.Callback = @CancelFcn;
    hcancel.Units = 'normalized';
    hcancel.Position = [0.8 0.03 0.15 0.1];

% OK
    hok = uicontrol('Style','pushbutton');
    hok.BackgroundColor = bgclr;
    hok.String = 'OK';
    hok.Callback = @OkFcn;
    hok.Units = 'normalized';
    hok.Position = [0.63 0.03 0.15 0.1];


    movegui(hwin,'center');
    set(hwin,'Visible','on');

    handles = guidata(hwin);
    handles.in_params = param_fit;
    handles.cancel_change = cancel_change;
    guidata(hwin,handles)

%% Callbacks
% add row
    function RowFcn(source,eventdata)
        handles = guidata(hwin);
        if isempty(handles.in_params)
            data = zeros(1,6);
        else
            data = handles.in_params;
            data(end+1,:) = data(end,:);
        end
        set(htable, 'data', data);
        handles.in_params = data;
        guidata(hwin,handles);
    end
% remove row
    function RowRemFcn(source,eventdata)
        handles = guidata(hwin);
        data = handles.in_params;
        data = data(1:end-1,:);
        set(htable, 'data', data);
        handles.in_params = data;
        guidata(hwin,handles);
    end
% OK
    function OkFcn(source,eventdata)
        handles = guidata(hwin);
        handles.in_params = htable.Data;
        guidata(hwin,handles);
        
        uiresume(hwin);
        close(hwin);
    end

%Cancel
    function CancelFcn(source,eventdata)
        handles = guidata(hwin);
        handles.in_params = handles.cancel_change;
        guidata(hwin,handles);
        
        uiresume(hwin);
        close(hwin);
    end

    uiwait(hwin)
    new_params = handles.in_params;
end