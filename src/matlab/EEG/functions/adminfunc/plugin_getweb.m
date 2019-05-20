% plugin_getweb - support function to get plugin information from the web

% Copyright (C) 2012- Arnaud Delorme
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function plugin = plugin_getweb(type, pluginOri, mode)

if nargin < 1, help plugin_getweb; return; end
if nargin < 2, pluginOri = []; end
if nargin < 3, mode = 'merge'; end % 'merge' or 'newlist' or 'update'

% convert plugin list format if necessary
if isfield(pluginOri, 'plugin'), pluginOri = plugin_convert(pluginOri); end

try
    disp( [ 'Retreiving plugins with extensions...' ] );
    if strcmpi(type, 'import')
    	[webPage, status] = plugin_urlread('https://sccn.ucsd.edu/wiki/Plugin_list_import');
    elseif strcmpi(type, 'process')
    	[webPage, status] = plugin_urlread('https://sccn.ucsd.edu/wiki/Plugin_list_process');
    else
    	[webPage, status] = plugin_urlread('https://sccn.ucsd.edu/wiki/Plugin_list_all');
    end
catch
    error('Cannot connect to the Internet to retrieve extension list');
end

% retreiving download statistics
try
    disp( [ 'Retreiving download statistics...' ] );
    [stats, status] = plugin_urlread('http://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_getcountall_with_rating.php');
    stats = textscan(stats, '%s%d%s%s%f%d', 'delimiter', char(9));
catch
    stats = {};
    disp('Cannot connect to the Internet to retrieve statistics for extensions');
end

if status == 0
    error('Cannot connect to the Internet to retrieve extension list');
end

% parse the web page
% ------------------
try
    plugin = parseTable(webPage);
catch
    error('Cannot parse extension list - please contact eeglab@sccn.ucsd.edu');
end

% find correspondance with plugin list
% ------------------------------------
if ~isempty(pluginOri)
     currentNames = lower({ pluginOri.name });
else currentNames = {};
end
allMatch = [];
allMatchPluginInd = [];
for iRow = 1:length(plugin)
    % fix links
    if isfield(plugin, 'zip'), plugin(iRow).zip = strrep(plugin(iRow).zip, '&amp;', '&'); end
        
    % get number of downloads
    if ~isempty(stats)
        indMatch = strmatch(lower(plugin(iRow).name), lower(stats{1}), 'exact');
        if ~isempty(indMatch)
             plugin(iRow).downloads = stats{2}(indMatch(1));
             if length(stats) > 2 && ~isempty(stats{4}{indMatch(1)})
                 plugin(iRow).version   = stats{3}{indMatch(1)};
                 plugin(iRow).zip       = stats{4}{indMatch(1)};
             end
             plugin(iRow).rating    = stats{5}(indMatch(1));
             plugin(iRow).numrating = stats{6}(indMatch(1));
        else
            plugin(iRow).downloads = 0;
        end
    else
        plugin(iRow).downloads = 0;
    end
    
    % match with existiting plugins
    indMatch = strmatch(lower(plugin(iRow).name), currentNames, 'exact');
    if isempty(indMatch)
        plugin(iRow).currentversion  = '-';
        plugin(iRow).installed       = 0;
        plugin(iRow).installorupdate = 1;
        plugin(iRow).status          = 'notinstalled';
    else
        if length(indMatch) > 1
            disp([ 'Warning: duplicate extension ' plugin(iRow).name ' instaled' ]); 
        end
        plugin(iRow).currentversion = pluginOri(indMatch).currentversion;
        plugin(iRow).foldername     = pluginOri(indMatch).foldername;
        plugin(iRow).status         = pluginOri(indMatch).status;
        plugin(iRow).installed      = 1;
        if strcmpi(plugin(iRow).currentversion, plugin(iRow).version)
            plugin(iRow).installorupdate = 0;
        else
            plugin(iRow).installorupdate = 1;
        end
        allMatch = [ allMatch indMatch(:)' ];
        allMatchPluginInd = [allMatchPluginInd iRow];
    end
    
end

if strcmpi(mode, 'update')
    plugin = plugin(allMatchPluginInd);
end

% put all the installed plugins first
% -----------------------------------
if ~isempty(plugin)
    [tmp reorder] = sort([plugin.installed], 'descend');
    plugin = plugin(reorder);
%     plugin(1).currentversion  = '0.9';
%     plugin(1).version         = '1';
%     plugin(1).foldername      = 'test';
%     plugin(1).installed       = 1;
%     plugin(1).installorupdate = 1;
%     plugin(1).description     = 'test';
%     plugin(1).webdoc          = 'test';
%     plugin(1).name            = 'test';
end

if strcmpi(mode, 'merge') && ~isempty(pluginOri)
    indices = setdiff([1:length(pluginOri)], allMatch);
    fields  = fieldnames(pluginOri);
    lenPlugin = length(plugin);
    
    for indPlugin = 1:length(indices)
        for indField = 1:length(fields)
            value  = getfield(pluginOri, { indices(indPlugin)  }, fields{ indField });
            plugin = setfield(plugin   , { lenPlugin+indPlugin }, fields{ indField }, value);
        end
    end
end

% parse the web table
% ===================
function plugin = parseTable(tmp);

plugin = [];
if isempty(tmp), return; end

% get table content
% -----------------
tableBeg = findstr('Plug-in name', tmp);
tableEnd = findstr('</table>', tmp(tableBeg:end));
tableContent = tmp(tableBeg:tableBeg+tableEnd-2);
endFirstLine = findstr('</tr>', tableContent);
tableContent = tableContent(endFirstLine(1)+5:end);

% parse table entries
% -------------------
posBegRow = findstr('<tr>' , tableContent);
posEndRow = findstr('</tr>', tableContent);
if length(posBegRow) ~= length(posEndRow) || isempty(posBegRow)
    error('Cannot connect to the Internet to retrieve plugin list');
end
for iRow = 1:length(posBegRow)
    rowContent = tableContent(posBegRow(iRow)+4:posEndRow(iRow)-1);
    posBegCol = findstr('<td>' , rowContent);
    posEndCol = findstr('</td>', rowContent);
    for iCol = 1:length(posBegCol)
        table{iRow,iCol} = rowContent(posBegCol(iCol)+4:posEndCol(iCol)-1);
    end
end

%% extract zip link and plugin name from first column
% --------------------------------------------------
%href="http://www.unicog.org/pm/uploads/MEG/ADJUST_PLUGIN.zip" class="external text" title="http://www.unicog.org/pm/uploads/MEG/ADJUST_PLUGIN.zip" rel="nofollow">ADJUST PLUGIN</a></td
for iRow = 1:size(table,1)
    
    % get link
    [plugin(iRow).name, plugin(iRow).webdoc] = parsehttplink(table{iRow,1});
    plugin(iRow).version = table{iRow,2};
    
    % description
    tmp = deblank(table{iRow,3}(end:-1:1));
    plugin(iRow).description = deblank(tmp(end:-1:1));
    
    % zip file
    [~,plugin(iRow).zip] = parsehttplink(table{iRow,4});

    % tags
    tmpTags = textscan(table{iRow,5}, '%s', 'delimiter', ',');
    tmpTags = tmpTags{1}';
    plugin(iRow).tags = tmpTags;
    tmpTags(2,:) = { ', ' };
    plugin(iRow).rawtags = [ tmpTags{:} ];
    plugin(iRow).rawtags(end-1:end) = [];
    
    % rating link
    plugin(iRow).webrating = [ 'https://sccn.ucsd.edu/eeglab/plugin_uploader/simplestar.php?plugin=' plugin(iRow).name '&version=' plugin(iRow).version ];
end

function [txt link] = parsehttplink(currentRow)
    openTag  = find(currentRow == '<');
    closeTag = find(currentRow == '>');
    if isempty(openTag)
        link = '';
        txt = currentRow;
    else
        % parse link
        link = currentRow(openTag(1)+1:closeTag(1)-1);
        hrefpos = findstr('href', link);
        link = link(hrefpos:end);
        quoteInd = find(link == '"');
        link = link(quoteInd(1)+1:quoteInd(2)-1);
        
        for iTag = length(openTag):-1:1
            currentRow(openTag(iTag):closeTag(iTag)) = [];
        end
        txt = currentRow;
    end

function plugin = plugin_convert(pluginOri)

for iRow = 1:length(pluginOri)
    plugin(iRow).currentversion = pluginOri(iRow).version;
    plugin(iRow).foldername     = pluginOri(iRow).foldername;
    plugin(iRow).status         = pluginOri(iRow).status;
    plugin(iRow).name           = pluginOri(iRow).plugin;
    plugin(iRow).installed      = 1;
end
