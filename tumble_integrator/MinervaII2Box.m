classdef MinervaII2Box < handle
    %MinervaII2 Geometry of MinervaII-2 rover for orientation plotting.
    %   This class is a dimension-accurate representation of the
    %   MinervaII-2 rover for orientation visualization purposes.
    %
    %   Inspired by the MATLAB HelperBox class
    %
    %   Due to the MATLAB copyright, this class should not be used for any
    %   official research or commercial applications.
    %
    %   Andrew Price
    %   21日01月2020年
    
    properties (Hidden)
        ColorIndex = 1;
    end
    
    properties
        AppWindow;
        
        pAxes;                  % Image axix' settings

        pRoverPatch;            % Patches (plotted panels of the rover)
        pRoverLines;            % Plotted line segments for the antennas
        pLineGroups;            % Antenna groupings
        pSolarPatch;            % Patches (plotted solar panels)
        pInitialPanelVertices;  % Initial locations of the panel vertices
        pInitialAntennaVertices;% Initial locations of the antenna vertices
        pInitialSolarVertices;  % Initial location of the solar panel vertices
        
        % figure default properties
        Title = '';
        AxesPosition = [0.1300 0.1100 0.7750 0.8150];
        XLimits = [-105, 105];
        YLimits = [-105, 105];
        ZLimits = [-105, 105];
        
        INITIAL_VIEW_ANGLE = [-1, 0, 0];
        
        pLinkedPropsObject
    end
    
    methods
        % Constructor
        function obj = MinervaII2Box(fig,varargin)
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
            %INITIALIZE Draw the MinervaII-2 rover geometry
            %   Initial orientation is
            %       x axis (panel 3) is point north
            %       y axis (panel 5) is pointing east
            %       z axis (panel 10) is pointing down
            
            % measurements from the MinervaII-2
            o1 = 57.4;      %[mm]
            o2 = 40.6;
            w = 2*o2+o1;    %[138.6 mm]
            %wa = 154.51;   % unused measured dimensions
            %wb = 172.55;
            h = 135.01;
            %a1 = 129.58;
            a2 = 4.5;
            a3 = 14;
            %ac = 236.96;
            %c = 14;
            %d = 150;
            %s1 = 22;        % solar panel dimensions
            %s2 = 22.5;
            %s3 = 43;
            %s4 = 33;
            
            % body
            xx = [0 , 0    , o2, o1+o2, w    , w , o1+o2, o2 ,0 , 0    , o2, o1+o2, w    , w  , o1+o2, o2]';
            yy = [o2, o1+o2, w , w    , o1+o2, o2, 0    , 0  ,o2, o1+o2, w , w    , o1+o2, o2 , 0    , 0 ]';
            zz = [0 , 0    , 0 , 0    , 0    , 0 , 0    , 0  ,h , h    , h , h    , h    , h  , h    , h ]';
            % antenna
            ax = [a2,a2,w-a2,w-a2,0.5*w,0.5*w,0.5*w,0.5*w,a2,a2,w-a2,w-a2,0.5*w,0.5*w,0.5*w,0.5*w]'; %vertice X coordinates
            ay = [0.5*w,0.5*w,0.5*w,0.5*w,a2,a2,w-a2,w-a2,0.5*w,0.5*w,0.5*w,0.5*w,a2,a2,w-a2,w-a2]'; %vertice Y coordinates
            az = [0,-a3,-a3,0,0,-a3,-a3,0,h,h+a3,h+a3,h,h,h+a3,h+a3,h]';                             %vertice Z coordinates
            % solar panels (geometry saved elsewhere, load it)
            Vars = {'GroupS', 'GSPV'};
            temp = load('SolarPanelPatchData.mat',Vars{:});
            obj.pInitialSolarVertices = temp.GSPV;
            
            % recentre the geometric centre to the origin
            offset = [0.5*w, 0.5*w, 0.5*h];
            xx = xx - offset(1);
            ax = ax - offset(1);
            yy = yy - offset(2);
            ay = ay - offset(2);
            zz = zz - offset(3);
            az = az - offset(3);
            
            % Construct the vertice matrices
            xxyyzz = [xx, yy, zz];
            axayaz = [ax, ay, az];
            obj.pInitialPanelVertices = xxyyzz;
            obj.pInitialAntennaVertices = axayaz;
                
            % Panel vertice vector indices
            panels  = [1,2,10,9,NaN,NaN,NaN,NaN;...     % panel 1
                2,3,11,10,NaN,NaN,NaN,NaN;...           % panel 2
                3,4,12,11,NaN,NaN,NaN,NaN;...           % panel 3
                4,5,13,12,NaN,NaN,NaN,NaN;...           % panel 4
                5,6,14,13,NaN,NaN,NaN,NaN;...           % panel 5
                6,7,15,14,NaN,NaN,NaN,NaN;...           % panel 6
                7,8,16,15,NaN,NaN,NaN,NaN;...           % panel 7
                8,1,9,16,NaN,NaN,NaN,NaN;...            % panel 8
                1,2,3,4,5,6,7,8;...                     % bottom
                9,10,11,12,13,14,15,16];                % top

            %          TOP DOWN VIEW
            %            Panel 3
            %             ------
            % Panel 2    /      \   Panel 4
            %          /          \
            % Panel 1 |     Top    | Panel 5
            %         |    Panel   |
            %          \          /
            % Panel 8    \      /   Panel 6
            %             ------
            %            Panel 7
            
            % Antenna vertice vector indices
            Groups = [1, 2, 3, 4;...     % bottom antenna 1
                5, 6, 7, 8;...         	% bottom antenna 2
                9, 10, 11, 12;...      	% top antenna 1
                13, 14, 15, 16];...    	% top antenna 2
            obj.pLineGroups = Groups;
            
            % Create the figure axis that all plots will use
            axx = createAxes(obj);
            
            % Colour data
            cdata = repmat(axx.ColorOrder(obj.ColorIndex,:), size(panels,1), 1);
            
            % Construct the panel geometry
            obj.pRoverPatch = patch(axx, 'Vertices', xxyyzz, 'Faces', panels,...
                'FaceVertexCData', cdata, 'FaceColor', 'flat');
            hold on     % Don't overwrite and delete the objects
            % Construct the antenna geometry
            L1 = plot3(axx,ax(Groups(1,:)),ay(Groups(1,:)),az(Groups(1,:)),...
                'Color','k','Marker','none','LineWidth',5, 'HandleVisibility', 'off');
            L2 = plot3(axx,ax(Groups(2,:)),ay(Groups(2,:)),az(Groups(2,:)),...
                'Color','k','Marker','none','LineWidth',5, 'HandleVisibility', 'off');
            L3 = plot3(axx,ax(Groups(3,:)),ay(Groups(3,:)),az(Groups(3,:)),...
                'Color','k','Marker','none','LineWidth',5, 'HandleVisibility', 'off');
            L4 = plot3(axx,ax(Groups(4,:)),ay(Groups(4,:)),az(Groups(4,:)),...
                'Color','k','Marker','none','LineWidth',5, 'HandleVisibility', 'off'); 
            obj.pRoverLines = [L1 L2 L3 L4];
            % Construct the solar panel geometry
            SP = cell(10,1);
            for i = 1:10
                SP{i} = patch(axx, 'Vertices', temp.GSPV{i}, 'Faces', temp.GroupS{i},...
                    'FaceVertexCData', [0 0 0], 'FaceColor', 'flat');
            end
            obj.pSolarPatch = [SP(1) SP(2) SP(3) SP(4) SP(5) SP(6) SP(7) SP(8) SP(9) SP(10)];
            
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
            RoverPanels = obj.pRoverPatch;
            Antenna = obj.pRoverLines;
            SolarPanels = obj.pSolarPatch;
            initPanelVertices = obj.pInitialPanelVertices;
            initAntennaVertices = obj.pInitialAntennaVertices;
            initSolarVertices = obj.pInitialSolarVertices;
            AGroups = obj.pLineGroups; % antenna vertice groups
            
            % Perform the rotation on the vertices
            newVertPanel =  initPanelVertices * R;
            newVertAntenna = initAntennaVertices * R;
            % newVertSolar = initSolarVertices * R;
            
            % Associate the new vertices to the objects
            set(RoverPanels, 'Vertices', newVertPanel);
            for i = 1:4
                set(Antenna(i),...
                    'XData',newVertAntenna(AGroups(i,:),1),...
                    'YData',newVertAntenna(AGroups(i,:),2),...
                    'ZData',newVertAntenna(AGroups(i,:),3));
            end
            for i = 1:10
                newVertSolar = initSolarVertices{i} * R;
                set(SolarPanels{i}, 'Vertices', newVertSolar);
            end
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

