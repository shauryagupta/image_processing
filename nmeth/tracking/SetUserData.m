function SetUserData (figHandle, object, replace, oName)
%SETUSERDATA connects any object to the figure 'figHandle'
%
% SYNOPSIS SetUserData(gHandle, object, oName)
%
% INPUT figHandle: a figure handle
%       object   : any object (number, string, struct...)
%       replace  : if replace = 1 the object is repalced in case of a name
%                  conflict. For any other value of replace an error occurs in
%                  case of a name conflict.
%       oName    : (optional) the name string under which the object is stored,
%                  by default the name of the 'object' variable is taken.
%
%This file is part of u-track.
%
%    u-track is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%    u-track is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with u-track.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright: 10/8/99	dT

uD = get(figHandle,'UserData');

% Check which name to use
if(nargin==4)
   name=oName;
else
   name=inputname(2);
end;

%Check if UserData is empty
if (isempty(uD))
   %initialize var
   uD=[];
else
   % Check for name conflict
   if ((any(strcmp(name,fieldnames(uD)))) & (replace~=1))
      error(['There exists already a field named ' 39 name 39]);
   end;
end;

uD=setfield(uD,name,object);
set(figHandle,'UserData',uD);
