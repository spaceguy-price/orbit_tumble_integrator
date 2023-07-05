classdef KHI_target < handle
    %KHI_target Geometry of KHI project for orientation plotting.
    %   Inspired by the MATLAB HelperBox class
    %
    %   Andrew Price
    %   10日01月2021年
    
    properties (Hidden)
        ColorIndex = 1;
    end
    
    properties
        AppWindow;
        
        pAxes;               % Image axix' settings

        pMainSidePatch;      % Patches (plotted panels of the rover)
        pMainEndPatch;
        pMainPanelVertices;  % Initial locations of the panel vertices
        pThrusterSidePatch; 
        pThrusterEndPatch;
        pThrusterPanelVertices;
        pNoseSidePatch;
        pNosePanelVertices;
        
        % figure default properties
        Title = '';
        AxesPosition = [0.1300 0.1100 0.7750 0.8150];
        XLimits = [-1, 1];
        YLimits = [-1, 1];
        ZLimits = [-1, 1];
        
        INITIAL_VIEW_ANGLE = [-0.01, -0.01, 0.01];
        
        pLinkedPropsObject
    end
    
    methods
        % Constructor
        function obj = KHI_target(fig,varargin)
            %MINERVAII2 Construct an instance of this class
            %   varargin: axesposition, title, colorindex
            obj.AppWindow = fig;
            if numel(varargin) > 0
                obj.AxesPosition = varargin{1};
            end
            if numel(varargin) > 1
                    obj.Title = varargin{2};
            end
            if numel(varargin) > 2
                obj.ColorIndex = varargin{3};
            end
            
            initialize(obj); %draw the initial object
        end
        
        function initialize(obj)
            %INITIALIZE 
            % Draw the KHI new target geometry
            % refer to MinervaIIBox.m for more object creation
            
            % -----1----- Construct the main cylinder
            rad = linspace(0,2*pi(),36);
            r = 0.25; %radius [m]
            Xc = zeros(2*length(rad),1);
            Yc = Xc;
            Zc = Xc;
            for i = 1:length(rad)
                Xc(i)               = r*cos(rad(i));
                Xc(i + length(rad)) = r*cos(rad(i));
                Yc(i)               = r*sin(rad(i));
                Yc(i + length(rad)) = r*sin(rad(i));
            end
            %Xc = Xc;                        %CoM x is located at 0 mm
            %Yc = Yc;                        %CoM y is located at 0 mm
            Zc(1:length(rad))     = - 0.5;   %CoM z is located at 500 mm
            Zc(length(rad)+1:end) =   0.5;
            
            % -----2----- Construct the thruster cylinder
            rad = linspace(0,2*pi(),36);
            r = 0.2; %radius [m]
            Xt = zeros(2*length(rad),1);
            Yt = Xt;
            Zt = Xt;
            for i = 1:length(rad)
                Xt(i)               = r*cos(rad(i));
                Xt(i + length(rad)) = r*cos(rad(i));
                Yt(i)               = r*sin(rad(i));
                Yt(i + length(rad)) = r*sin(rad(i));
            end
            %Xt = Xt;                        
            %Yt = Yt;                        
            Zt(1:length(rad))     = - 0.7;   
            Zt(length(rad)+1:end) = - 0.5;
            
            % -----3----- Construct the nose cone
            rad = linspace(0,2*pi(),36);
            r = 0.25; %radius [m]
            Xn = zeros(2*length(rad),1);
            Yn = Xn;
            Zn = Xn;
            for i = 1:length(rad)
                Xn(i)               = r*cos(rad(i));
                Xn(i + length(rad)) = 0;
                Yn(i)               = r*sin(rad(i));
                Yn(i + length(rad)) = 0;
            end
            %Xn = Xn;                       
            %Yn = Yn;                       
            Zn(1:length(rad))     = 0.5;    
            Zn(length(rad)+1:end) = 0.8;
            
            % Construct the vertice matrices
            xxyyzzc = [Xc, Yc, Zc];
            xxyyzzt = [Xt, Yt, Zt];
            xxyyzzn = [Xn, Yn, Zn];
            obj.pMainPanelVertices = xxyyzzc;
            obj.pThrusterPanelVertices = xxyyzzt;
            obj.pNosePanelVertices = xxyyzzn;
                
            % Panel vertice vector indices
            temp1 = [1:1:length(rad)]';
            temp2 = temp1 + 1;
            temp2(end) = 1;
            temp4 = [length(rad)+1:1:2*length(rad)]';
            temp3 = temp4 + 1;
            temp3(end) = length(rad) + 1;
            
            panels = [temp1, temp2, temp3, temp4];
            clear temp1 temp2 temp3 temp4 f
            
            ends = [1:1:length(rad); length(rad)+1:1:2*length(rad)];
            
            % Create the figure axis that all plots will use
            axx = createAxes(obj);
            
            % Colour data
            cdata = repmat(axx.ColorOrder(obj.ColorIndex,:), size(panels,1), 1);
            cdata2 = cdata;
            cdata(1,:) = [1,0,0];   %set the first bar to red
            cdata(10,:)= [0,1,0];   %set the 90deg bar to green
            cdata(19,:)= [0.5,0,0]; %set the 180 bar light red
            cdata(28,:)= [0,0.75,0]; %set the 270 bar light green
            cdata3 = repmat(axx.ColorOrder(obj.ColorIndex,:), size(ends,1), 1);
            
            % Construct the panel geometry
            obj.pMainSidePatch = patch(axx, 'Vertices', xxyyzzc, 'Faces',...
                panels,'FaceVertexCData', cdata, 'FaceColor', 'flat');
            hold on     % Don't overwrite and delete the objects
            obj.pMainEndPatch = patch(axx, 'Vertices', xxyyzzc, 'Faces',...
                ends,'FaceVertexCData', cdata3, 'FaceColor', 'flat');
            
            obj.pThrusterSidePatch = patch(axx, 'Vertices', xxyyzzt, 'Faces',...
                panels,'FaceVertexCData', cdata2, 'FaceColor', 'flat');
            obj.pThrusterEndPatch = patch(axx, 'Vertices', xxyyzzt, 'Faces',...
                ends,'FaceVertexCData', cdata3, 'FaceColor', 'flat');
            
            obj.pNoseSidePatch = patch(axx, 'Vertices', xxyyzzn, 'Faces',...
                panels,'FaceVertexCData', cdata, 'FaceColor', 'flat');
        end
        
        function axx = createAxes(obj)
            %CREATEAXES Draw basic axes for box plotting

            fig = obj.AppWindow;

            axx = axes(fig, 'OuterPosition', obj.AxesPosition);
            
            axx.Title.String = obj.Title;

            view(axx,obj.INITIAL_VIEW_ANGLE);  

            % axis equal;
            axx.DataAspectRatioMode = 'manual';
            axx.DataAspectRatio = [1 1 1];
            axx.PlotBoxAspectRatioMode = 'manual';
            axx.PlotBoxAspectRatio = [1.2 1 1];
            
            % Reference frame is NED. Reverted all directions to normal to
            % make the visualization process easier.
            axx.XDir = 'normal';            
            axx.YDir = 'normal';    %'reverse'
            axx.ZDir = 'normal';    %'reverse'
            axx.XLabel.String = 'x (North)';
            axx.YLabel.String = 'y (East)';
            axx.ZLabel.String = 'z (Down)';

            axx.XGrid = 'on';
            axx.XMinorGrid = 'on';
            axx.YGrid = 'on';
            axx.YMinorGrid = 'on';
            axx.ZGrid = 'on';
            axx.ZMinorGrid = 'on';

            axx.XLimMode = 'manual';
            axx.XLim = obj.XLimits;
            axx.YLimMode = 'manual';
            axx.YLim = obj.YLimits;
            axx.ZLimMode = 'manual';
            axx.ZLim = obj.ZLimits;
            
            obj.pAxes = axx;
        end
        
        function deleteAxes(obj)
            delete([obj.pAxes]);
        end
        
        function update(obj, R)
            %UPDATE - update the orientation of rover geometry
            if isa(R, 'quaternion')         % check if the rotation information is a quaternion
                R = rotmat(R, 'frame');     % if it is a quaternion, convert to a rotation matrix
            end
            
            % Pull the geometry data from the properties
            mainSidePanels          = obj.pMainSidePatch;
            mainEndPanels           = obj.pMainEndPatch;
            mainPanelVertices       = obj.pMainPanelVertices;
            thrusterSidePanels      = obj.pThrusterSidePatch;
            thrusterEndPanels       = obj.pThrusterEndPatch;
            thrusterPanelVertices   = obj.pThrusterPanelVertices;
            noseSidePanels          = obj.pNoseSidePatch;
            nosePanelVertices       = obj.pNosePanelVertices;
            
            % Perform the rotation on the vertices
            newMainPanel     = mainPanelVertices * R;
            newThrusterPanel = thrusterPanelVertices * R;
            newNosePanel     = nosePanelVertices * R;
            
            % Associate the new vertices to the objects
            set(mainSidePanels,     'Vertices', newMainPanel);
            set(mainEndPanels,      'Vertices', newMainPanel);
            set(thrusterSidePanels, 'Vertices', newThrusterPanel);
            set(thrusterEndPanels,  'Vertices', newThrusterPanel);
            set(noseSidePanels,     'Vertices', newNosePanel);
        end
        
        function syncRotationCallback(objs)
        %SYNCROTATIONCALLBACK - link rotation to another set of axes
            linkedPropsObject = linkprop([objs.pAxes], {'XLim', 'YLim', 'ZLim', ...
                'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'});
            for i = 1:numel(objs)
                objs(i).pLinkedPropsObject = linkedPropsObject;
            end
        end
    end
end

