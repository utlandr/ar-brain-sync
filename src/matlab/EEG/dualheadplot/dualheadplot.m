% dualheadplot() - plot a spherically-splined EEG field map on a semi-realistic
%                                 3-D head model. Plot the links of a
%                                 connectivity matrix
%
% Standard-mode Usage thereafter:
%       >> [hdaxis cbaraxis] =
%       dualheadplot('spline_file','Param','Value',...)
%       Example:
%       nChans = 64
%       A = triu(rand(nChans),1);
%       A = A+A';     % this is a random connectivity matrix (values between [0 1] as synchrony/coherence)
%                            %NOTICE that the diagonal is formed by zeros
%       var2BePloted = sum(A);   % values to be plotted at each electrode position (e.g. the number of connections, the power, etc)
%      figure,
%       %to plot the variables (eg. the power) on the scalp, without conexion links
%      dualheadplot('aurelieDataKarimHead.spl','values2plot',var2BePloted,'maplimits',[min(var2BePloted),max(var2BePloted)],'electrodes','on','cbar',0,'labels',2);
%
%       %To plot the variables (eg. the power) on the scalp, with the conexion links ()eg synchrony, coherence, etc whose values are between 0.99 and 1
%       dualheadplot('aurelieDataKarimHead.spl','adjacencyMatrix',A,'thresholdLink',[0.99 1],'values2plot',var2BePloted,'maplimits',[min(var2BePloted),max(var2BePloted)],'electrodes','on','cbar',0,'labels',0);
%
%
% Required Standard-mode Inputs:
%
%   'spline_file' - spline filename, computed and saved in 'setup' mode (see eeglab and pop_headplot)
%
% Optional Standard-mode Inputs:
%
%   'values2plot'  - vector containing a data value at each electrode position
%                       {default zeros(1,numberOfElectrodes)}
%   'adjacencyMatrix'  - the connectivity (e.g. synchrony/coherence) matrix
%   
%                                IMPORTANT: there are only zeros on the diagonal
%                               {default  zeros(numberOfElectrodes,numberOfElectrodes)}
%   ''thresholdLink''  - real valued vector with the threshold for plotting the links {default [1.01 1.1]}
%   'meshfile'   - [string] mesh file name. See file content in the setup-mode
%                  description above. {default: the EEGLAB head template file}.
%   'electrodes' - ['on'|'off'] -> show electrode positions {default 'on'}
%   'title'      -  Plot title {default: none}
%   'labels'     -  2 -> plot stored electrode labels;
%                   1 -> plot channel numbers; 0 -> no labels {default 0}
%   'cbar'       -  0 -> Plot colorbar {default: no colorbar}
%                        Note: standard jet colormap) red = +;blue = -;green=0
%                   h -> Colorbar axis handle (to specify headplot location)
%   'view'       - Camera viewpoint in deg. [azimuth elevation]
%                  'back'|'b'=[  0 30]; 'front'|'f'=[180 30]
%                  'left'|'l'=[-90 30]; 'right'|'r'=[ 90 30];
%                  'frontleft'|'bl','backright'|'br', etc.,
%                  'top'=[0 90],  Can rotate with mouse {default [143 18]}
%   'maplimits'  - 'absmax' -> make limits +/- the absolute-max
%                  'maxmin' -> scale to data range
%                   [min,max] -> user-definined values
%                      {default = 'absmax'}
%   'lights'     - (3,N) matrix whose rows give [x y z] pos. of each of
%                   N lights {default: four lights at corners}
%   'lighting'   - 'off' = show wire frame head {default 'on'}
%   'colormap'   -  3-column colormap matrix {default: jet(64)}
%   'verbose'    - 'off' -> no msgs, no rotate3d {default: 'on'}
%
% Outputs:
%      hdaxis    - handle of the head axes (the existing gca when headplot() was called)
%      cbaraxis  - handle of the color bar axes
%
% Note: if an error is generated, dualheadplot() may close the current figure
%
% ---------------------------------------------------------------------------------------
% Original Authors: Arnaud Delorme, Colin Humphries, Scott Makeig, SCCN/INC/UCSD,
%                           La Jolla, 1998-
%   Modified from headplot()  EEGLAB
%
% ---------------------------------------------------------------------------------------
% Plotarcs routine originally produced by: Miguel Valencia, CNRS / LENA, 2007
% Embedded and simplified by: Mario Chavez, CNRS/LENA, 2010
% Hyperscanning version by: Guillaume Dumas, CNRS/LENA 2012 & Institut Pasteur 2016
% ---------------------------------------------------------------------------------------
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [HeadAxes, ColorbarHandle] = dualheadplot(arg1, varargin)
%[HeadAxes, ColorbarHandle] = dualheadplot(values, AdjMtrx, thr, arg1, varargin)
if nargin < 1
    help dualheadplot
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Set Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFAULT_MESH = ['mheadnew.mat'];      % nice head for the Lille data (size 501 KB)

DEFAULT_LIGHTS = [-125  125  80; ...
    125  125  80; ...
    125 -125 125; ...
    -125 -125 125; ...% ];    % default lights at four corners
      0     10   -80]; % add light from below!

HeadCenter = [0 0 0];
FaceColor  = [.8 .55 .35]*1.2; % ~= ruddy Caucasian  --->0 more mexican color !!!

plotelecopt.NamesDFac  = 1.0;  % plot electrode names/numbers out from markers
plotelecopt.NamesColor = 'k';
plotelecopt.NamesSize  =  12;   % FontSize for electrode names
plotelecopt.MarkerSize = 20;
plotelecopt.MarkerColor= 'k';

sqaxis     = 1;     % if non-zero, make head proportions anatomical
title_font = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin <2
    disp('Warning: Only the electrodes will be displayed !  Type help dualheadplot')
end
spline_file = arg1;

g = finputcheck( varargin, { ...
    'values2plot'   'real'    []        [];  %vector, one value per electrode (default = zeros)
    'adjacencyMatrix'  'real'  []    [];  %the connectivity matrix (default = zeroes matrix)
    'thresholdLink'     'real'   []    [1.01   1.1]; % threshold for plotting thelinks
    'cbar'       'real'   [0 Inf]         []; % Colorbar value must be 0 or axis handle.
    'lighting'   'string' { 'on' 'off' }  'on';
    'verbose'    'string' { 'on' 'off' }  'on';
    'maplimits'  { 'string' 'real' }  []  'absmax';
    'title'      'string' []              '';
    'lights'     'real'   []              DEFAULT_LIGHTS;
    'view'       'real'   []              [143 18];
    'colormap'   'real'   []              jet(64);
    'transform'  'real'   []              [];
    'meshfile'   'string' []              DEFAULT_MESH;
    'electrodes' 'string' { 'on' 'off' }  'on';
    'orilocs'    { 'string' 'struct' } [] '';
    'labels'     'integer' [0 1 2]        0 }, 'headplot');
if isstr(g) error(g); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open head mesh and electrode spline files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(spline_file)
    error(sprintf('headplot(): spline_file "%s" not found. Run headplot in "setup" mode\n',...
        spline_file));
end
load(spline_file, '-mat');

if ~isempty(g.values2plot)
    values = g.values2plot;
else
    values = zeros(1,length(Xe)*2);
end

%indices=indices-2; %FLAG
if exist('indices'),
    try,
        values = values([indices indices+length(values)/2]);
    catch, error('problem of index or electrode number with splinefile'); end
end
enum = length(values);
if enum ~= length(Xe)*2
    close;
   error(sprintf('%s: Number of values in spline file should equal number of electrodes\n',mfilename));
end


% load mesh file
% --------------
if ~exist(g.meshfile)    
    error(sprintf('%s: mesh file "%s" not found\n',mfilename, g.meshfile));
end
load(g.meshfile,'-mat');

if exist('index1') ~= 1, index1 = sort(unique(TRI1(:))); end
if exist('TRI2')   ~= 1, TRI2 = []; end
if exist('NORM')   ~= 1, NORM = []; end
if exist('TRI1')   ~= 1, error('Variable ''TRI1'' not defined in mesh file'); end
if exist('POS')    ~= 1, error('Variable ''POS'' not defined in mesh file'); end
if exist('center') ~= 1, center = [0 0 0]; disp('Using [0 0 0] for center of head mesh'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%

meanval = mean(values); values = values - meanval; % make mean zero
onemat = ones(enum,1);
lamd = 0.1;

C1 = pinv([(G + lamd);ones(1,enum/2)]) * [values(1:length(values)/2)';0]; % fixing division error
P1 = zeros(1,size(gx,1));
for j = 1:size(gx,1)
    P1(j) = dot(C1,gx(j,:));
end
P1 = P1 + meanval;

C2 = pinv([(G + lamd);ones(1,enum/2)]) * [values(length(values)/2+1:length(values))';0]; % fixing division error
P2 = zeros(1,size(gx,1));
for j = 1:size(gx,1)
    P2(j) = dot(C2,gx(j,:));
end
P2 = P2 + meanval;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%
cla % clear axis
HeadAxes = gca;
P=[P1 P2];
W = zeros(1,size(POS,1));
m = size(g.colormap,1);
if size(g.maplimits) == [1,2]
    amin = g.maplimits(1);
    amax = g.maplimits(2);
elseif strcmp(g.maplimits,'maxmin') | strcmp(g.maplimits,'minmax')
    amin = min(min(abs(P)))*1.02; % 2% shrinkage keeps within color bounds
    amax = max(max(abs(P)))*1.02;
elseif strcmp(g.maplimits,'absmax')
    amin = min(min(abs(P)))*1.02; % 2% shrinkage keeps within color bounds
    amax = max(max(abs(P)))*1.02;
    amax = max(-amin, amax);
    amin = -amax;
    if amax == amin
        amax = amin+eps;
    end
    %amin = -max(max(abs(P)))*1.02; % 2% shrinkage keeps within color bounds
    %amax = -amin;
end

idx1 = min(m,round((m-1)*(P(1:length(P)/2)-amin)/(amax-amin))+1); % get colormap indices
idx2 = min(m,round((m-1)*(P(length(P)/2+1:length(P))-amin)/(amax-amin))+1); % get colormap indices
%idx = round((m-1)*P/(amax-amin))+m/2;
%idx = max(1,min(m,idx)); % get colormap indices

W1(index1) = idx1;
W2(index1) = idx2;
colormap(g.colormap)
POS=cat(2,POS(:,2),-POS(:,1),POS(:,3));
%p1 = patch('Vertices',POS,'Faces',TRI1,'FaceVertexCdata',W(:),...
%    'FaceColor','flat', 'cdatamapping', 'scaled', 'tag', 'mesh');    %%%%%%%%% Plot scalp map %%%%%%%%%
p1 = patch('Vertices',POS+cat(2,ones(size(POS,1),1)*75,zeros(size(POS,1),1),zeros(size(POS,1),1)),'Faces',TRI1,'FaceVertexCdata',W1(:),...
    'FaceColor','interp', 'cdatamapping', 'direct', 'tag', 'mesh');    %%%%%%%%% Plot scalp map %%%%%%%%%

POS2=cat(2,POS(:,1)*cos(pi)+POS(:,2)*sin(pi)+300,POS(:,1)*sin(pi)+POS(:,2)*cos(pi),POS(:,3));
p2 = patch('Vertices',POS2+cat(2,ones(size(POS2,1),1)*75,zeros(size(POS2,1),1),zeros(size(POS2,1),1)),'Faces',TRI1,'FaceVertexCdata',W2(:),...
    'FaceColor','interp', 'cdatamapping', 'direct', 'tag', 'mesh');    %%%%%%%%% Plot scalp map %%%%%%%%%

if exist('NORM') == 1 & ~isempty(NORM)
   set(p1, 'vertexnormals', NORM);
   set(p2, 'vertexnormals', NORM);
end

if ~isempty(TRI2)
    FCmap = [g.colormap; g.colormap(end,:); FaceColor; FaceColor; FaceColor];
    colormap(FCmap)
    W = ones(1,size(POS,1))*(m+2);
    p3 = patch('Vertices',POS+cat(2,ones(size(POS,1),1)*75,zeros(size(POS,1),1),zeros(size(POS,1),1)),'Faces',TRI2,'FaceColor','interp',...
        'FaceVertexCdata',W(:)); %%%%%%%% Plot face and lower head %%%%%%
    POS2=cat(2,POS(:,1)*cos(pi)+POS(:,2)*sin(pi)+300,POS(:,1)*sin(pi)+POS(:,2)*cos(pi),POS(:,3));
    p4 = patch('Vertices',POS2+cat(2,ones(size(POS2,1),1)*75,zeros(size(POS2,1),1),zeros(size(POS2,1),1)),'Faces',TRI2,'FaceColor','interp',...
        'FaceVertexCdata',W(:)); %%%%%%%% Plot face and lower head %%%%%%
else
    p3 = [];
    p4 = [];
end

%axis([-125 125 -125 125 -125 125])
axis off % hide axis frame

%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw colorbar - Note: uses enhanced cbar() function by Colin Humphries
%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(g.cbar)
    BACKCOLOR = get(gcf,'Color');
    if g.cbar == 0
        ColorbarHandle = MYcbar(0,3,[amin amax]);
        pos = get(ColorbarHandle,'position');  % move left & shrink to match head size
        set(ColorbarHandle,'position',[pos(1)-.05 pos(2)+0.13 pos(3)*0.7 pos(4)-0.26]);
    else
        ColorbarHandle = MYcbar(g.cbar,3,[amin amax]);
    end
end
axes(HeadAxes);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Turn on lights
%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(g.lighting,'on')
    set([p1 p2 p3 p4],'EdgeColor','none')
    
    for i = 1:size(g.lights,1)
        hl(i) = light('Position',g.lights(i,:),'Color',0.5*[1 1 1],...
            'Style','infinite');
    end
    if ~isempty(p3)
        set(p2,'DiffuseStrength',.6,'SpecularStrength',1,...
            'AmbientStrength',.4,'SpecularExponent',15,'SpecularColorReflectance',0.3,'BAckFAceLighting', 'lit')
    end
    if ~isempty(p4)
        set(p2,'DiffuseStrength',.6,'SpecularStrength',1,...
            'AmbientStrength',.4,'SpecularExponent',15,'SpecularColorReflectance',0.3,'BAckFAceLighting', 'lit')
    end
    set(p1,'DiffuseStrength',.6,'SpecularStrength',1,...
        'AmbientStrength',.3,'SpecularExponent',15,'SpecularColorReflectance',0.3,'BAckFAceLighting', 'lit')
    set(p2,'DiffuseStrength',.6,'SpecularStrength',1,...
        'AmbientStrength',.3,'SpecularExponent',15,'SpecularColorReflectance',0.3,'BAckFAceLighting', 'lit')
    lighting phong  % all this gives a matte reflectance
    %lighting gouraud 
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Set viewpoint
%%%%%%%%%%%%%%%%%%%%%%%%%

if isstr(g.view)
    switch lower(g.view)
        case {'front','f'}
            view(-180,30)
        case {'back','b'}
            view(0,30)
        case {'left','l'}
            view(-90,30)
        case {'right','r'}
            view(90,30)
        case {'frontright','fr'}
            view(135,30)
        case {'backright','br'}
            view(45,30)
        case {'frontleft','fl'}
            view(-135,30)
        case {'backleft','bl'}
            view(-45,30)
        case 'top'
            view(0,90)
        case 'bottom'    % undocumented option!
            view(0,-90)
            Lights = [-125 125 80;   ...
                125 125 80;   ...
                125 -125 125; ...
                -125 -125 125; ...
                0 10 -80]; % add light from below!
        otherwise
            close; error(sprintf('%s: Invalid View value %s\n',mfilename,g.view));
    end
else
    if ~isstr(g.view)
        [h,a] = size(g.view);
        if h~= 1 | a~=2
            close; error(sprintf('%s: View matrix size must be (1,2).',mfilename));
        end
    end
    %camlight headlight
    view(g.view)   % set camera viewpoint
end

thr = g.thresholdLink;
if ~isempty(g.adjacencyMatrix)
    AdjMtrx = g.adjacencyMatrix;
else
    AdjMtrx = zeros(size(newElect,1),size(newElect,1));
end

if exist('center') 
    HeadCenter= center; 
end
HeadCenter=mean(POS);
newElect(:,1)=-newElect(:,1);

plotarcs(AdjMtrx, thr, newElect, HeadCenter, plotelecopt); % plot the arcs locations

if strcmp(g.electrodes,'on') % plot the electrode locations
    if exist('newElect')
        plotelecopt.labelflag = g.labels;
        plotelec(newElect, ElectrodeNames, HeadCenter, plotelecopt);
    else
        fprintf('Variable newElect not read from spline file.\n');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turn on rotate3d, allowing rotation of the plot using the mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(g.verbose,'on')
    rotate3d on;   % Allow 3-D rotation of the plot by dragging the
else             % left mouse button while cursor is on the plot
    rotate3d off
end
% Make axis square
if sqaxis
    axis image    % keep the head proportions human and as large as possible
end
% Add a plot title
if ~isempty(g.title);
    % title(['\n' g.title],'fontsize',title_font);
    title([g.title],'fontsize',title_font); % Note: \n not interpreted by matlab-5.2
end
%cameramenu
return

% %%%%%%%%%%%%%%%
% plot electrodes
% %%%%%%%%%%%%%%%
function plotelec(newElect, ElectrodeNames, HeadCenter, opt);
newElect=cat(2,newElect(:,2),newElect(:,1),newElect(:,3));
newNames1 = newElect*opt.NamesDFac; % Calculate electrode label positions
newNames2 = cat(2,newNames1(:,1)*cos(pi)+newNames1(:,2)*sin(pi)+300,newNames1(:,1)*sin(pi)+newNames1(:,2)*cos(pi),newNames1(:,3));
newNames1 = newNames1+cat(2,ones(size(newNames1,1),1)*75,zeros(size(newNames1,1),1),zeros(size(newNames1,1),1));
newNames2 = newNames2+cat(2,ones(size(newNames2,1),1)*75,zeros(size(newNames2,1),1),zeros(size(newNames2,1),1));

fprintf('%s: Interpolation for %i  electrodes: \n', mfilename,length(newElect));

for i = 1:size(newElect,1)
    if newElect(i,:) ~= [0 0 0]  % plot radial lines to electrode sites
        %         line([newElect(i,1) HeadCenter(1)],[newElect(i,2) HeadCenter(2)],...
        %             [newElect(i,3) HeadCenter(3)],'color',opt.MarkerColor,'linewidth',1);
        
        if opt.labelflag == 1        % plot electrode numbers
            t1=text(newNames1(i,1),newNames1(i,2),newNames1(i,3),int2str(i));
            set(t1,'Color',opt.NamesColor,'FontSize',opt.NamesSize,'FontWeight','bold',...
                'HorizontalAlignment','center');
            t2=text(newNames2(i,1),newNames2(i,2),newNames2(i,3),int2str(i));
            set(t2,'Color',opt.NamesColor,'FontSize',opt.NamesSize,'FontWeight','bold',...
                'HorizontalAlignment','center');
        elseif opt.labelflag == 2   % plot electrode names
            if exist('ElectrodeNames')
                name = sprintf('%s',ElectrodeNames(i,:));
                t1=text(newNames1(i,1),newNames1(i,2),newNames1(i,3),name);
                set(t1,'Color',opt.NamesColor,'FontSize',opt.NamesSize,'FontWeight','bold',...
                    'HorizontalAlignment','center');
                t2=text(newNames2(i,1),newNames2(i,2),newNames2(i,3),name);
                set(t2,'Color',opt.NamesColor,'FontSize',opt.NamesSize,'FontWeight','bold',...
                    'HorizontalAlignment','center');
            else
                fprintf('Variable ElectrodeNames not read from spline file.\n');
            end
        else               % plot electrode markers
            line(newNames1(i,1),newNames1(i,2),newNames1(i,3),'marker',...
                '.','markersize',opt.MarkerSize,'color',opt.MarkerColor,'linestyle','none');
            line(newNames2(i,1),newNames2(i,2),newNames2(i,3),'marker',...
                '.','markersize',opt.MarkerSize,'color',opt.MarkerColor,'linestyle','none');
        end
    end
end


% %%%%%%%%%%%%%%%%%%%%%
% plot relation arcs
% %%%%%%%%%%%%%%%%%%%%%
function plotarcs(AdjMtrx, thr, newElect, HeadCenter, opt);
%newElect = newElect*opt.NamesDFac;
newElect=cat(2,newElect(:,2),newElect(:,1),newElect(:,3));
newNames1 = newElect*opt.NamesDFac; % Calculate electrode label positions
newNames2 = cat(2,newNames1(:,1)*cos(pi)+newNames1(:,2)*sin(pi)+300,newNames1(:,1)*sin(pi)+newNames1(:,2)*cos(pi),newNames1(:,3));
newNames1 = newNames1+cat(2,ones(size(newNames1,1),1)*75,zeros(size(newNames1,1),1),zeros(size(newNames1,1),1));
newNames2 = newNames2+cat(2,ones(size(newNames2,1),1)*75,zeros(size(newNames2,1),1),zeros(size(newNames2,1),1));
newElect=cat(1,newNames1,newNames2);
hold on

thr1 = min(thr);
thr2 = max(thr);

[I,J] = find(((AdjMtrx>thr1).*(AdjMtrx<thr2))> 0);

cmap = jet(128);cmap=cmap(65:128,:);
nColors = size(cmap,1);

for n = 1:length(I)
    if (I(n)<=size(AdjMtrx,1)/2)&&(J(n)<=size(AdjMtrx,1)/2)
        originArcs=mean(newNames1)-[0 0 60];
    else
        if (I(n)>size(AdjMtrx,1)/2)&&(J(n)>size(AdjMtrx,1)/2)
            originArcs=mean(newNames2)-[0 0 60];
        else
            originArcs=mean(newElect)-[0 0 60];
        end
    end
    arc = plotarc([newElect(I(n),1),newElect(I(n),2),newElect(I(n),3)],...
    [newElect(J(n),1),newElect(J(n),2),newElect(J(n),3)],originArcs);%[185.5 0 -75.5]);
    %        [newElect(J(n),1),newElect(J(n),2),newElect(J(n),3)],[0 0 0]);

    value = AdjMtrx(I(n),J(n));
    idx = round(nColors*(value-thr1)/(thr2-thr1)); if idx < 1, idx = 1; end
    h = plot3(arc(1,:),arc(2,:),arc(3,:));
%     set(h,'LineWidth',10*(value-thr1)/(thr2-thr1),'Color',cmap(idx,:))
        set(h,'LineWidth',5,'Color',cmap(idx,:)) %change by alej after commenting previous
end

[I,J] = find(((AdjMtrx<-thr1).*(AdjMtrx>-thr2))> 0);

cmap = jet(128);cmap=cmap(1:64,:);
nColors = size(cmap,1);

for n = 1:length(I)
    if (I(n)<=size(AdjMtrx,1)/2)&&(J(n)<=size(AdjMtrx,1)/2)
        originArcs=mean(newNames1)-[0 0 60];
    else
        if (I(n)>size(AdjMtrx,1)/2)&&(J(n)>size(AdjMtrx,1)/2)
            originArcs=mean(newNames2)-[0 0 60];
        else
            originArcs=mean(newElect)-[0 0 60];
        end
    end
    arc = plotarc([newElect(I(n),1),newElect(I(n),2),newElect(I(n),3)],...
    [newElect(J(n),1),newElect(J(n),2),newElect(J(n),3)],originArcs);%[185.5 0 -75.5]);%
    %        [newElect(J(n),1),newElect(J(n),2),newElect(J(n),3)],[0 0 0]);

    value = AdjMtrx(I(n),J(n));
    idx = round(nColors*(-value-thr1)/(thr2-thr1)); if idx < 1, idx = 1; end
    h = plot3(arc(1,:),arc(2,:),arc(3,:));
%     set(h,'LineWidth',10*(-value-thr1)/(thr2-thr1),'Color',cmap(idx,:))
        set(h,'LineWidth',5,'Color',cmap(idx,:)) %change by alej after commenting previous
end

% %%%%%%%%%%%%%%%%%%%%%
% Estimate the arcs
% %%%%%%%%%%%%%%%%%%%%%
function arc = plotarc(p1,p2,c);

if nargin < 1
    p1 = [-4 6 2];
    p2 = [4 5 2];
    c = [0 0 10];
end

pm = mean([p1;p2]);

v1 = [p2(1)-p1(1), p2(2)-p1(2), p2(3)-p1(3)];
v1 = v1/norm(v1);

v2 = [pm(1)-c(1), pm(2)-c(2), pm(3)-c(3)];
v2 = v2/norm(v2);

v3 = cross(v2,v1);
v3 = v3/norm(v3);

u1 = v1/norm(v1);

u2 = v2 - dot(v2,u1) * u1;
u2 = u2/norm(u2);

u3 = v3 - dot(v3,u1) * u1 - dot(v3,u2) * u2;
u3 = u3/norm(u3);

A = [[1,0,0];[0,1,0];[0,0,1]];
B = [u1;u2;u3];

C = inv(B')*A';

radius = sqrt((p2(1)-p1(1))^2 + (p2(2)-p1(2))^2 + (p2(3)-p1(3))^2)/2;
angl = 0.0:4:180;
d2r = pi/180;
for i = 1:length(angl)
    Arc_(1,i) = radius*cos(angl(i)*d2r);
    Arc_(2,i) = radius*sin(angl(i)*d2r);
    Arc_(3,i) = 0;
end

arc = C' * Arc_;
arc(1,:) = arc(1,:)+pm(1);
arc(2,:) = arc(2,:)+pm(2);
arc(3,:) = arc(3,:)+pm(3);


function [g, varargnew] = finputcheck( vararg, fieldlist, callfunc, mode )
% finputcheck() - check Matlab function {'key','value'} input argument pairs
%
% Usage: >> result = finputcheck( varargin, fieldlist );
%        >> [result varargin] = finputcheck( varargin, fieldlist, ... 
%                                                         callingfunc, mode );
% Input:
%   varargin  - Cell array 'varargin' argument from a function call using 'key', 
%               'value' argument pairs. See Matlab function 'varargin'
%   fieldlist - A 4-column cell array, one row per 'key'. The first
%               column contains the key string, the second its type(s), 
%               the third the accepted value range, and the fourth the 
%               default value.  Allowed types are 'boolean', 'integer', 
%               'real', 'string', 'cell' or 'struct'.  For example,
%                       {'key1' 'string' { 'string1' 'string2' } 'defaultval_key1'}
%                       {'key2' {'real' 'integer'} { minint maxint } 'defaultval_key2'} 
%  callingfunc - Calling function name for error messages. {default: none}.
%  mode        - ['ignore'|'error'] ignore keywords that are either not specified 
%                in the fieldlist cell array or generate an error. 
%                {default: 'error'}.
% Outputs:
%   result     - If no error, structure with 'key' as fields and 'value' as 
%                content. If error this output contain the string error.
%   varargin   - residual varagin containing unrecognized input arguments.
%                Requires mode 'ignore' above.
%
% Note: In case of error, a string is returned containing the error message
%       instead of a structure.
%
% Example (insert the following at the beginning of your function):
%	result = finputcheck(varargin, ...
%               { 'title'         'string'   []       ''; ...
%                 'percent'       'real'     [0 1]    1 ; ...
%                 'elecamp'       'integer'  [1:10]   [] });
%   if isstr(result)
%       error(result);
%   end
%
% Note: 
%   The 'title' argument should be a string. {no default value}
%   The 'percent' argument should be a real number between 0 and 1. {default: 1}
%   The 'elecamp' argument should be an integer between 1 and 10 (inclusive).
%
%   Now 'g.title' will contain the title arg (if any, else the default ''), etc.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 July 2002


% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 10 July 2002, arno@salk.edu

	if nargin < 2
		help finputcheck;
		return;
	end
	if nargin < 3
		callfunc = '';
	else 
		callfunc = [callfunc ' ' ];
	end
    if nargin < 4
        mode = 'do not ignore';
    end
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
	
	varargnew = {};
	% create structure
	% ----------------
	if ~isempty(vararg)
		for index=1:length(vararg)
			if iscell(vararg{index})
				vararg{index} = {vararg{index}};
			end
		end
		try
			g = struct(vararg{:});
		catch
            vararg = removedup(vararg);
            try,
                g = struct(vararg{:});
            catch
                g = [ callfunc 'error: bad ''key'', ''val'' sequence' ]; return;
            end
		end
	else 
		g = [];
	end
	
	for index = 1:size(fieldlist,NAME)
		% check if present
		% ----------------
		if ~isfield(g, fieldlist{index, NAME})
			g = setfield( g, fieldlist{index, NAME}, fieldlist{index, DEF});
		end
		tmpval = getfield( g, {1}, fieldlist{index, NAME});
		
		% check type
		% ----------
        if ~iscell( fieldlist{index, TYPE} )
            res = fieldtest( fieldlist{index, NAME},  fieldlist{index, TYPE}, ...
                           fieldlist{index, VALS}, tmpval, callfunc );
            if isstr(res), g = res; return; end
        else 
            testres = 0;
            tmplist = fieldlist;
            for it = 1:length( fieldlist{index, TYPE} )
                if ~iscell(fieldlist{index, VALS})
                     res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}, tmpval, callfunc );
                else res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}{it}, tmpval, callfunc );
                end
                if ~isstr(res{it}), testres = 1; end
            end
            if testres == 0,
                g = res{1};
                for tmpi = 2:length(res)
                    g = [ g 10 'or ' res{tmpi} ];
                end
                return; 
            end
        end
	end
    
    % check if fields are defined
	% ---------------------------
	allfields = fieldnames(g);
	for index=1:length(allfields)
		if isempty(strmatch(allfields{index}, fieldlist(:, 1)', 'exact'))
			if ~strcmpi(mode, 'ignore')
				g = [ callfunc 'error: undefined argument ''' allfields{index} '''']; return;
			end
			varargnew{end+1} = allfields{index};
			varargnew{end+1} = getfield(g, {1}, allfields{index});
		end
	end


function g = fieldtest( fieldname, fieldtype, fieldval, tmpval, callfunc );
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
    g = [];
    
    switch fieldtype
     case { 'integer' 'real' 'boolean' 'float' }, 
      if ~isnumeric(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be numeric' ]; return;
      end
      if strcmpi(fieldtype, 'boolean')
          if tmpval ~=0 & tmpval ~= 1
              g = [ callfunc 'error: argument ''' fieldname ''' must be 0 or 1' ]; return;
          end  
      else 
          if strcmpi(fieldtype, 'integer')
              if ~isempty(fieldval)
                  if (isnan(tmpval) & ~any(isnan(fieldval))) ...
                          & (~ismember(tmpval, fieldval))
                      g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
                  end
              end
          else % real or float
              if ~isempty(fieldval)
                  if tmpval < fieldval(1) | tmpval > fieldval(2)
                      g = [ callfunc 'error: value out of range for argument ''' fieldname '''' ]; return;
                  end
              end
          end
      end  
      
      
     case 'string'
      if ~isstr(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a string' ]; return;
      end
      if ~isempty(fieldval)
          if isempty(strmatch(lower(tmpval), lower(fieldval), 'exact'))
              g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
          end
      end

      
     case 'cell'
      if ~iscell(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a cell array' ]; return;
      end
      
      
     case 'struct'
      if ~isstruct(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a structure' ]; return;
      end
      
      
     case '';
     otherwise, error([ 'finputcheck error: unrecognized type ''' fieldname '''' ]);
    end

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
% make sure if all the values passed to unique() are strings, if not, exist
%try
    [tmp indices] = unique(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        fprintf('Note: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end
    cella = cella(sort(union(indices*2-1, indices*2)));
%catch
    % some elements of cella were not string
%    error('some ''key'' values are not string.');
%end    
return


function [handle]=MYcbar(arg,colors,minmax, grad)
% cbar() - Display full or partial color bar
%
% Usage:
%    >> cbar % create a vertical cbar on the right side of a figure
%    >> cbar(type) % specify direction as 'vert' or 'horiz'
%    >> cbar(type,colors) % specify which colormap colors to plot
%  else
%    >> cbar(axhandle) % specify the axes to draw cbar in
%
%    >> h = cbar(type|axhandle,colors, minmax, grad)
%
% Inputs:
%  type      - ['vert'|'horiz'] direction of the cbar {default: 'vert')
%              ELSE axhandle = handle of axes to draw the cbar
%  colors    - vector of colormap indices to display, or integer to truncate upper
%              limit by.
%              (int n -> display colors [1:end-n]) {default: 0}
%  minmax    - [min, max] range of values to label on colorbar
%  grad      - [integer] number of tick labels. {default: 5}.
%
% Example:
%         >> colormap('default') % default colormap is 64-color 'jet'
%         >> cbar('vert',33:64); % plot a vertical cbar colored green->red
%                                % useful for showing >0 (warm) and 0 (green)
%                                % values only in a green=0 plot
%
% Author: Colin Humphries, Arnaud Delorme, CNL / Salk Institute, Feb. 1998-
%
% See also: colorbar()



if nargin < 2
    colors = 0;
end
posscale = 'off';
if nargin < 1
    arg = 'vert';
    ax = [];
else
    if isempty(arg)
        arg = 0;
    end
    if arg(1) == 0
        ax = [];
        arg = 'vert';
    elseif strcmpi(arg, 'pos')
        ax = [];
        arg = 'vert';
        posscale = 'on';
    else
        if ischar(arg)
            ax = [];
        else
            ax = arg;
            arg = [];
        end
    end
end

if nargin>2
    if size(minmax,1) ~= 1 | size(minmax,2) ~= 2
        help cbar
        fprintf('cbar() : minmax arg must be [min,max]\n');
        return
    end
end
if nargin < 4
    grad = 5;
end

%obj = findobj('tag','cbar','parent',gcf);
%if ~isempty(obj) & ~isempty(arg)
%  arg = [];
%  ax = obj;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose colorbar position
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (length(colors) == 1) & (colors == 0)
    t = caxis;
else
    t = [0 1];
end
if ~isempty(arg)
    if strcmp(arg,'vert')
        cax = gca;
        pos = get(cax,'Position');
        stripe = 0.04;
        edge = 0.01;
        space = .02;
        
        %    set(cax,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
        %    rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
        
        set(cax,'Position',[pos(1) pos(2) pos(3) pos(4)])
        rect = [pos(1)+pos(3)+space pos(2) stripe*pos(3) pos(4)];
        ax = axes('Position', rect);
    elseif strcmp(arg,'horiz')
        cax = gca;
        pos = get(cax,'Position');
        stripe = 0.075;
        space = .1;
        set(cax,'Position',...
            [pos(1) pos(2)+(stripe+space)*pos(4) pos(3) (1-stripe-space)*pos(4)])
        rect = [pos(1) pos(2) pos(3) stripe*pos(4)];
        ax = axes('Position', rect);
    end
else
    pos = get(ax,'Position');
    if pos(3) > pos(4)
        arg = 'horiz';
    else
        arg = 'vert';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw colorbar using image()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map = colormap;
n = size(map,1);

if length(colors) == 1
    if strcmp(arg,'vert')
        if strcmpi(posscale, 'on')
            image([0 1],[0 t(2)],[ceil(n/2):n-colors]');
        else
            image([0 1],t,[1:n-colors]');
        end
        set(ax,'xticklabelmode','manual')
        set(ax,'xticklabel',[],'YAxisLocation','right')
        
    else
        image(t,[0 1],[1:n-colors]);
        set(ax,'yticklabelmode','manual')
        set(ax,'yticklabel',[],'YAxisLocation','right')
    end
    set(ax,'Ydir','normal','YAxisLocation','right')
    
else % length > 1
    
    if max(colors) > n
        error('Color vector excedes size of colormap')
    end
    if strcmp(arg,'vert')
        image([0 1],t,[colors]');
        set(ax,'xticklabelmode','manual')
        set(ax,'xticklabel',[])
    else
        image([0 1],t,[colors]);
        set(ax,'yticklabelmode','manual')
        set(ax,'yticklabel',[],'YAxisLocation','right')
    end
    set(ax,'Ydir','normal','YAxisLocation','right')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust cbar ticklabels
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 2
    if strcmp(arg,'vert')
        Cax = get(ax,'Ylim');
    else
        Cax = get(ax,'Xlim');
    end
    CBTicks = [Cax(1):(Cax(2)-Cax(1))/(grad-1):Cax(2)]; % caxis tick positions
    CBLabels = [minmax(1):(minmax(2)-minmax(1))/(grad-1):minmax(2)]; % tick labels
    
    dec = floor(log10(max(abs(minmax)))); % decade of largest abs value
    CBLabels = ([minmax]* [ linspace(1,0, grad);linspace(0, 1, grad)]);
    %[1.0 .75 .50 .25 0.0; 0.0 .25 .50 .75 1.0]);
    if dec<1
        CBLabels = round(CBLabels*10^(1-dec))*10^(dec-1);
    elseif dec == 1
        CBLabels = round(CBLabels*10^(2-dec))*10^(dec-2);
    else
        CBLabels = round(CBLabels);
    end
    % minmax
    % CBTicks
    % CBLabels
    
    if strcmp(arg,'vert')
        set(ax,'Ytick',CBTicks);
        set(ax,'Yticklabel',CBLabels);
    else
        set(ax,'Xtick',CBTicks);
        set(ax,'Xticklabel',CBLabels);
    end
end
handle = ax;

%%%%%%%%%%%%%%%%%%
% Adjust cbar tag
%%%%%%%%%%%%%%%%%%

set(ax,'tag','cbar')
return