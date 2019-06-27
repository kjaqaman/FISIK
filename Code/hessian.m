function varargout=hessian(data)
% HESSIAN compute the hessian of data
%
% SYNOPSIS  [FXX, FXY,FXZ,FYX,FYY,FYZ,FZX,FZY,FZZ]=curvature3D(img,pnt)
%
% INPUT data   : 2D or 3D data
% OUTPUT  FXX...FZZ   : hessian entries for all positions
%
% Copyright (C) 2019, Jaqaman Lab - UTSouthwestern
%
% This file is part of FISIK.
% 
% FISIK is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% FISIK is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with FISIK.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% c: 7/12/01	dT
delta=1;
if (ndims(data)==3)
    [FX,FY,FZ] = gradient(data,delta);
    [FXX,FXY,FXZ] = gradient(FX,delta);
    [FYX,FYY,FYZ] = gradient(FY,delta);
    [FZX,FZY,FZZ] = gradient(FZ,delta);
    varargout(1:9)={FXX,FXY,FXZ,FYX,FYY,FYZ,FZX,FZY,FZZ};
elseif(ndims(data)==3)
    [FX,FY] = gradient(data,delta);
    [FXX,FXY] = gradient(FX,delta);
    [FYX,FYY] = gradient(FY,delta);
    varargout(1:4)={FXX,FXY,FYX,FYY};
end;
