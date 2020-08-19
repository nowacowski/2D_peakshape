function [fit_results, fit_errors, tw] = Fit2D_gui(data)
%{
    GUI program for peakshape analysis using multivariate normal distribution.
    Program takes stack of two-dimensional arrays, with its 3 axis and fits
    sum of specified number of 2D Gaussian functions:

    f(x,y) = A * exp(- 1 / (2 * (1 - rho^2)) * ((x - mu_x)^2 / sigma_x^2 + (y - mu_y)^2 / sigma_y^2 - 2*rho*(x - mu_x)*(y - mu_y) / sigma_x*sigma_y)
    
    where A - amplitude, mu_x/mu_y - peak center in x/y axis,
    sigma_x/sigma_y - width of the peak in x/y axis, rho - correlation
    between x and y axes

    Then saves fitting parameters with errors as csv file.

    Recommended usage:
    - specify slices which will be fitted between 'Tw start' and 'Tw stop'
    - choose region of interest for X and Y axes 'Get ROI'
        - click on graph and drag cursor for required ROI
        - confirm by double click on edge of drawn ROI
    - specify offset between 0:1
        - offset 0.2 means that data below 0.2*max(abs(data)) will not be
        included in fitting procedure
    - fit to data
        - for single 2D Gaussian click 'Start fitting'
            - program will perform fitting analytically without the need of
            giving initial parameters
        - for sum of 2 or more Gaussians specify initial parameters
        (performed iteratively)
            - click 'Initial parameters' and fill the table
            - each new row is an additional Gaussian function added
    - to inspect fitting look at the graphs in the middle panel
        - switch between 2D data slices using 'Tw' dropdown
        - use sliders to look at slices in each axis of specified 2D data

    * if you think fitting is not satisfactonary you can change initial
    parameters and run again, table will be filled with fitting results
    from the fitting to first data slice
    * provided initial parameters are used for all fitted slices
    * saved data for each 2D slice is a row of form: 
        [Tw, A_1, error, mu_x_1, error, sigma_x_1, error, mu_y_1, error,
        sigma_y_1, error, rho_1, error, ..., A_n, error, mu_x_n, error, sigma_x_n, error, mu_y_n, error,
        sigma_y_n, error, rho_n, error]
%}
%%
if isstruct(data) ~= 1
    error('Input data needs to be a structure array')
end
% check if axis are sorted in ascending order
if issorted(data.X, 'ascend') ~= 1
    error('X needs to be sorted in ascending order')
end
if issorted(data.Y, 'ascend') ~= 1
    error('Y needs to be sorted in ascending order')
end

X = data.X;
Y = data.Y;
T = data.T;
data = data.Data;

%%
% Initialization

Tw = floor(length(T)/2);
T_roi = T;

xROI = 1:length(X);
yROI = 1:length(Y);

xx = X(xROI(floor(length(xROI)/2)));
yy = Y(yROI(floor(length(yROI)/2)));

offset = 0;

param_fit = [];

path = [];
in_param = [];

% Main window
hwin = figure('Visible','off');
hwin.NumberTitle = 'off';
hwin.Name = '2D Gaussian fit';
hwin.Resize = 'on';
hwin.ResizeFcn = @WinResizeFcn;
hwin.Units = 'normalized';
hwin.Position = [0.1 0.1 0.75 0.6];

bgclr = get(hwin,'Color');

% panel 1
uip1 = uipanel('Title', 'Data');
uip1.Units = 'normalized';
uip1.Position = [0 0 0.444 1];

% panel 2
uip2 = uipanel('Title', 'Fitting');
uip2.Units = 'normalized';
uip2.Position = [0.444 0 0.48 1];

% panel 3
uip3 = uipanel('Title', 'Save results');
uip3.Units = 'normalized';
uip3.Position = [0.9240 0 0.076 1];

%handles = guidata(hwin);
handles.fig = hwin;
handles.uipanel3 = uip1;
handles.uipanel3 = uip2;
handles.uipanel3 = uip3;
handles.xROI = xROI;
handles.yROI = yROI;
handles.path = path;
handles.in_param = in_param;
handles.param_fit = param_fit;
guidata(hwin,handles)
%% buttons
% Close
hclose = uicontrol('Style','pushbutton');
hclose.BackgroundColor = bgclr;
hclose.String = 'Close';
hclose.Callback = @CloseFcn;
hclose.Units = 'normalized';

% Save
hsave = uicontrol('Style','pushbutton');
hsave.BackgroundColor = bgclr;
hsave.String = 'Save';
hsave.Callback = @SaveFcn;
hsave.Units = 'normalized';

% get ROI
hgetROI = uicontrol('Style','pushbutton');
hgetROI.BackgroundColor = bgclr;
hgetROI.String = 'Get ROI';
hgetROI.Callback = @getROIFcn;
hgetROI.Units = 'normalized';
hgetROI.TooltipString = sprintf('Select ROI on graph.\n Double click on edge.');

% xROI
hxROI = uicontrol('Style','Text');
hxROI.BackgroundColor = bgclr;
hxROI.String = [num2str(round(X(xROI(1)),2)),' : ',num2str(round(X(xROI(end)),2))];
hxROI.Units = 'normalized';

hxROItxt = uicontrol('Style','text');
hxROItxt.String = 'xROI:';
hxROItxt.HorizontalAlignment = 'left';
hxROItxt.BackgroundColor = bgclr;
hxROItxt.Units = 'normalized';

% yROI
hyROI = uicontrol('Style','Text');
hyROI.BackgroundColor = bgclr;
hyROI.String = [num2str(round(Y(yROI(1)),2)),' : ',num2str(round(Y(yROI(end)),2))];
hyROI.Units = 'normalized';

hyROItxt = uicontrol('Style','text');
hyROItxt.String = 'yROI:';
hyROItxt.HorizontalAlignment = 'left';
hyROItxt.BackgroundColor = bgclr;
hyROItxt.Units = 'normalized';

% cell array of Tw strings
nT = length(T);
Tstr = cell(1,nT);
for iT = 1:nT
    Tstr{iT} = num2str(T(iT));
end

% Tw
hTw = uicontrol('Style','popupmenu');
hTw.String = Tstr;
hTw.Callback = @TwFcn;
hTw.Units = 'normalized';
hTw.Value = length(T)/2;

hTwtxt = uicontrol('Style','text');
hTwtxt.String = 'Tw:';
hTwtxt.HorizontalAlignment = 'left';
hTwtxt.BackgroundColor = bgclr;
hTwtxt.Units = 'normalized';

% Tw start
hTstart = uicontrol('Style','popupmenu');
hTstart.String = Tstr;
hTstart.Units = 'normalized';
hTstart.Value = 1;

hTstarttxt = uicontrol('Style','text');
hTstarttxt.String = 'Tw start:';
hTstarttxt.HorizontalAlignment = 'left';
hTstarttxt.BackgroundColor = bgclr;
hTstarttxt.Units = 'normalized';

% Tw stop
hTstop = uicontrol('Style','popupmenu');
hTstop.String = Tstr;
hTstop.Value = length(T);
hTstop.Units = 'normalized';

hTstoptxt = uicontrol('Style','text');
hTstoptxt.String = 'Tw stop:';
hTstoptxt.HorizontalAlignment = 'left';
hTstoptxt.BackgroundColor = bgclr;
hTstoptxt.Units = 'normalized';

% offset
hoffset = uicontrol('Style','edit');
hoffset.Callback = @offsetFcn;
hoffset.Units = 'normalized';
hoffset.String = num2str(offset);

hofftxt = uicontrol('Style','text');
hofftxt.String = 'Offset:';
hofftxt.HorizontalAlignment = 'left';
hofftxt.BackgroundColor = bgclr;
hofftxt.Units = 'normalized';

% fitting
hfit = uicontrol('Style','pushbutton');
hfit.BackgroundColor = bgclr;
hfit.String = 'Start fitting';
hfit.Callback = @Fitting;
hfit.Units = 'normalized';

%change fitting parameters
hinparam = uicontrol('Style','pushbutton');
hinparam.BackgroundColor = bgclr;
hinparam.String = 'Initial Parameters';
hinparam.Callback = @InParam;
hinparam.Units = 'normalized';

% fitted y slider
hyAx = uicontrol('Style','slider');
hyAx.SliderStep = [0.01 0.01];
hyAx.Min = Y(yROI(1));
hyAx.Max = Y(yROI(end));
hyAx.Value = yy;
hyAx.SliderStep = [0.05 0.05];
hyAx.Callback = {@yAxFcn};
hyAx.Units = 'normalized';

% fitted x slider
hxAx = uicontrol('Style','slider');
hxAx.SliderStep = [0.01 0.01];
hxAx.Min = X(xROI(1));
hxAx.Max = X(xROI(end));
hxAx.Value = xx;
hxAx.SliderStep = [0.05 0.05];
hxAx.Callback = {@xAxFcn};
hxAx.Units = 'normalized';

% Move the GUI to the center of the screen
movegui(hwin,'center');

% Make the GUI visible
set(hwin,'Visible','on');

%% %%%%%% Callback functions %%%%%%%%%%%%

    fig1(X,Y,data,Tw);
    fig2(X,xROI,Y,yROI,data,xx,yy,param_fit)

%% Tw
    function TwFcn(source,eventdata)
        Tw = get(hTw,'Value');
        fig1(X,Y,data,Tw);
        fig2(X,xROI,Y,yROI,data,xx,yy,param_fit)        
    end

%% offset
    function offsetFcn(source,eventdata)
        offset = str2double(get(hoffset,'String'));
    end

%% getROI
    function getROIFcn(source,eventdata)
        fig1(X,Y,data,Tw);
        h = imrect;
        position = wait(h)
        
        [~,xROI1] = min((X-position(1)).^2);
        [~,xROI2] = min((X-(position(1)+position(3))).^2);
        xROI = xROI1:xROI2;
        
        [~,yROI1] = min((Y-position(2)).^2);
        [~,yROI2] = min((Y-(position(2)+position(4))).^2);
        yROI = yROI1:yROI2;
        
        set(hxROI,'String',[num2str(round(X(xROI(1)),2)),' : ',num2str(round(X(xROI(end)),2))]);
        set(hyROI,'String',[num2str(round(Y(yROI(1)),2)),' : ',num2str(round(Y(yROI(end)),2))]);
        
        xx = X(xROI(floor(length(xROI)/2)));
        yy = Y(yROI(floor(length(yROI)/2)));
        
        set(hyAx, 'min', Y(yROI(1)));
        set(hyAx, 'max', Y(yROI(end)));
        set(hxAx, 'min', X(xROI(1)));
        set(hxAx, 'max', X(xROI(end)));
        
        set(hxAx, 'Value', X(xROI(floor(length(xROI)/2))));
        set(hyAx, 'Value', Y(yROI(floor(length(yROI)/2))));

        fig2(X,xROI,Y,yROI,data,xx,yy,param_fit)
                
        handles = guidata(hwin);
        handles.xROI = xROI;
        handles.yROI = yROI;
        guidata(hwin,handles);
        
    end

%% Figure 1
    function fig1(X,Y,data,Tw)
   %             x   y   wx   wy
        pos = [0.25 0.1 0.73 0.74];
        
        axes(uip1);
        subplot('Position',pos(1,:));
        contourf(X,Y,data(:,:,Tw),20)
        shading interp     
    end

%% sliders for checking fit
    function yAxFcn(source,eventdata)
        yy = get(hyAx,'Value');
        fig2(X,xROI,Y,yROI,data,xx,yy,param_fit)
    end
    function xAxFcn(source,eventdata)
        xx = get(hxAx,'Value');
        fig2(X,xROI,Y,yROI,data,xx,yy,param_fit)
    end

%% Figure 2 / fit
    function fig2(X,xROI,Y,yROI,data,xx,yy,param_fit)
        
        Tcheck = get(hTw,'Value');
        Tstart = get(hTstart,'Value');
        
        [X2, Y2] = meshgrid(X(xROI),Y(yROI));
        xdata(:,:,1) = X2; 
        xdata(:,:,2) = Y2;
        border = ones(size(data(yROI,xROI,1)));
        
        % choose slice and Tw
        [~,nx] = min((xx-X(xROI)).^2);
        [~,ny] = min((yy-Y(yROI)).^2);
        
        % get fit values for graphs
        if ~isempty(param_fit) && (Tcheck-Tstart+1 >= 1)
            x = param_fit(Tcheck-Tstart+1,:);
            fit = Gauss2D(x,xdata,border);
            xslicefit = fit(ny,:);
            yslicefit = fit(:,nx);
        end
        
        % data for plots of slices
        xslice = data(yROI(1)+ny,xROI,Tcheck);
        yslice = data(yROI,xROI(1)+nx,Tcheck);
               
        %position of subplota
    %             x     y       wx      wy
        pos = [0.1    0.18    0.64     0.62;...
               0.1     0.81   0.64    0.15;...
               0.75    0.18    0.15     0.62];
        
        axes(uip2);
        % contourplot
        subplot('Position',pos(1,:));
        pcolor(X(xROI),Y(yROI),data(yROI,xROI,Tcheck));
        set(gca,'yticklabel',{[]},'xticklabel',{[]})
        shading interp
        if ~isempty(param_fit) && (Tcheck-Tstart+1 >= 1)
            hold on
            contour(X(xROI),Y(yROI),fit,16,'LineColor',[.8,.8,.8])
        end
        
        % horizontal line across contourplot
        hold on
        plot(X(xROI),interp1([X(xROI(1)),X(xROI(end))],[Y(yROI(ny)),Y(yROI(ny))],X(xROI)),'k')
        
        % vertical line across contourplot
        hold on
        plot(interp1([Y(yROI(end)),Y(yROI(1))],[X(xROI(nx)),X(xROI(nx))],Y(yROI)),Y(yROI),'k')
        hold off
        
        % horizontal slice plot
        subplot('Position',pos(2,:));
        if ~isempty(param_fit) && (Tcheck-Tstart+1 >= 1)
            plot(X(xROI),xslice,'o',X(xROI),xslicefit,'r')
        else
            plot(X(xROI),xslice,'o')
        end
        set(gca,'yticklabel',{[]},'XAxisLocation', 'top');%,'xticklabel',{[]})%,'ydir', 'reverse')
        xlim([X(xROI(1)) X(xROI(end))]);
        ylim([min(min(data(:,:,Tcheck))) max(max(data(yROI,xROI,Tcheck)))]);
        
        % vertical slice plot
        subplot('Position',pos(3,:));
        if ~isempty(param_fit) && (Tcheck-Tstart+1 >= 1)
            plot(yslice,Y(yROI),'o',yslicefit,Y(yROI),'r')
        else
            scatter(yslice,Y(yROI),'o')
        end
        set(gca,'xticklabel',{[]},'YAxisLocation', 'right');%,'yticklabel',{[]})%,'xdir', 'reverse')
        ylim([Y(yROI(1)) Y(yROI(end))]);
        xlim([min(min(data(:,:,Tcheck))) max(max(data(yROI,xROI,Tcheck)))]);
    end

%% Fitting
    function Fitting(source,eventdata)
        
        param_fit = [];
        Tstart = get(hTstart,'Value');
        Tstop = get(hTstop,'Value');
        T_roi = T(Tstart:Tstop);
        
        handles = guidata(hwin);
        xROI = handles.xROI;
        yROI = handles.yROI;
        offset = str2double(get(hoffset,'String'));
        
        % load provided initial parameters or estimate if there is no (only
        % for fitting one peak)
        if ~isempty(handles.in_param)
            x0 = handles.in_param;
        else
            % estimate initial parameters
            data_1 = data(yROI,xROI,Tstart);
            
            [~,mux0,wx0] = caruana(X(xROI),sum(data_1,1));  % fit gaussian analyticaly
            [~,muy0,wy0] = caruana(Y(yROI),sum(data_1,2));  % fit gaussian analyticaly
            data_1_v = data_1(:);
            [~, idx] = max(abs(data_1_v));
            A0 = data_1(idx);
            
            x0 = [A0, mux0, wx0, muy0, wy0, 0];
        end
        
        % boundries
        lb = zeros(size(x0));
        ub = zeros(size(x0));
        for i = 1:size(x0,1)
            if x0(i,1) > 0
                lb(i,:) = [x0(i,1:end-1)-x0*0.5, -1];
                ub(i,:) = [x0(i,1:end-1)+x0*0.5, 1];
            else
                lb(i,:) = [x0(i,1)+x0(i,1)*0.5, x0(i,2:end-1)-x0(i,2:end-1)*0.5, -1];
                ub(i,:) = [x0(i,1)-x0(i,1)*0.5, x0(i,2:end-1)+x0(i,2:end-1)*0.5, 1];
            end
        end
        % reshape initial parameters and boundries into vectors
        x0 = reshape(x0',1,[]);
        lb = reshape(lb',1,[]);
        ub = reshape(ub',1,[]);        
        
        % create meshgrid for axis
        [X2, Y2] = meshgrid(X(xROI),Y(yROI));
        xdata(:,:,1) = X2; 
        xdata(:,:,2) = Y2;
        
        %initialize progress bar
        bar = waitbar(1/length(T_roi),'Fitting progress'); 

        % fitting
        fit_err = zeros(length(T_roi), length(x0));
        for i = 1:length(T_roi)    %iterate for each T
            data_roi = data(yROI,xROI,Tstart+i-1);  %data slice
            border=abs(data_roi)>=(offset*(max(max(abs(data_roi)))));   %offset
            
            
            Options = optimset('TolFun',1e-9);
            
            fval = 1;
            old_f = 0;
            x = x0;
            while fval~=old_f      %repeat until the fit stop getting better
                old_f = fval;
                x0 = x;
                [x,fval,~,~,~,~,J] = lsqcurvefit(@(x,xdata)Gauss2D(x,xdata,border),x0,xdata,data_roi.*border,lb,ub,Options);
            end
            
            cov_err = cov_error(fval,xdata,x,J);
            fit_err(i,:) = sqrt(diag(cov_err));
            
            waitbar(i/length(T_roi),bar,'Fitting progress');
            param_fit(i,:) = x;
        end
        close(bar)  %stop progress bar
        
        handles = guidata(hwin);
        handles.param_fit = param_fit;
        handles.fit_err = fit_err; 
        handles.T_roi = T_roi;
        guidata(hwin,handles);
        fig2(X,xROI,Y,yROI,data,xx,yy,param_fit)
    end

%% Set initial parameters
    function InParam(source,eventdata)
        
        handles = guidata(hwin);
        if isempty(handles.param_fit)
            param_fit = handles.in_param;
        else
            param_fit = handles.param_fit(1,:);
        end
        
        in_param = InParamGui(param_fit);
        
        handles.in_param = in_param;
        guidata(hwin,handles);
    end

%%
% Close GUI
    function CloseFcn(source,eventdata)
        uiresume(gcbf);
        close(hwin);
    end
%%
% Save data
    function SaveFcn(source,eventdata)
        
        Tstart = get(hTstart,'Value');
        Tstop = get(hTstop,'Value');
        T_roi = T(Tstart:Tstop);
        handles = guidata(hwin);
        result = handles.param_fit;
        errors = handles.fit_err;
        
        dane = zeros(size(result,1),size(result,2)*2+1);
        dane(:,1) = T_roi;
        dane(:,2:2:end-1) = result;
        dane(:,3:2:end) = errors;
        
        [FileName,Path] = uiputfile({'*.csv','CSV Files (*.csv)';'*.*','All Files (*.*)'},'Save',path);
        path_file = fullfile(Path,FileName);
        
        [~, d] = version;
        if d > datetime(2019,1,1)
            writematrix(dane, path_file, 'Delimiter', ',');
        else
            csvwrite(path_file, dane);
        end
        
        handles = guidata(hwin);
        handles.path = Path;
        guidata(hwin,handles);
    end
%% Positions and sizes of elements
    function WinResizeFcn(source,eventdata)

        set(hclose,'Position',[0.94,0.45,0.05,0.08]);   %close button
        set(hsave,'Position',[0.94,0.6,0.05,0.08]);    %save button
        
        set(hTwtxt,'Position',[0.017,0.7,0.048,0.027]);   %Tw text
        set(hTw,'Position',[0.017,0.66,0.048,0.018]);      %Tw dropdown
        set(hTstarttxt,'Position',[0.017,0.52,0.048,0.027]);   %Tstart text
        set(hTstart,'Position',[0.017,0.48,0.048,0.018]);  %Tstart dropdown
        set(hTstoptxt,'Position',[0.017,0.34,0.048,0.027]);%Tstop text
        set(hTstop,'Position',[0.017,0.3,0.048,0.018]);   %Tstop dropdown

        set(hgetROI,'Position',[0.14,0.86,0.048,0.05]);    %getROI button
        set(hyROItxt,'Position',[0.2,0.9,0.048,0.03]);  %yROI text
        set(hyROI,'Position',[0.215,0.86,0.072,0.04]);     %yROI value
        set(hxROItxt,'Position',[0.3,0.9,0.048,0.03]);  %xROI text
        set(hxROI,'Position',[0.315,0.86,0.072,0.04]);     %xROI value
             
        set(hfit,'Position',[0.73,0.02,0.072,0.07]);    %fit button
        set(hinparam,'Position',[0.6,0.03,0.096,0.05]);     %initial parameters button
        set(hyAx,'Position',[0.47,0.15,0.018,0.65]);   %vertical slider
        set(hxAx,'Position',[0.479,0.11,0.3312,0.04]);   %horizontal slider

        set(hoffset,'Position',[0.53,0.03,0.036,0.05]);    %offset input
        set(hofftxt,'Position',[0.493,0.04,0.036,0.027]);%offset text
               
    end


uiwait(gcf);

% function outputs

fit_results = handles.param_fit;
fit_errors = handles.fit_err;
tw = handles.T_roi;
end