function ecef = ned2ecef(ned, orgllh)
% ecef2ned: converts from NED coordinates to ECEF coordinates.
%   writing by Liu yx, based on NaveGo
%   Copyright (C) 2014, Rodrigo Gonzalez, all rights reserved. 
%   
%   This file is part of NaveGo, an open-source MATLAB toolbox for 
%   simulation of integrated navigation systems.
%     
%   NaveGo is free software: you can redistribute it and/or modify
%   it under the terms of the GNU Lesser General Public License (LGPL) 
%   version 3 as published by the Free Software Foundation.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
% 
%   You should have received a copy of the GNU Lesser General Public 
%   License along with this program. If not, see 
%   <http://www.gnu.org/licenses/>.
%
% Reference: 
%			R. Gonzalez, J. Giribet, and H. Patiño. NaveGo: a 
% simulation framework for low-cost integrated navigation systems, 
% Journal of Control Engineering and Applied Informatics}, vol. 17, 
% issue 2, pp. 110-120, 2015. Inverse process of Eq. 15.
%
% Version: 001
% Date:    2016/0/11
% Author:  liuyx<liuyx@bupt.edu.cn>
% URL:     https://github.com/rodralez/navego 

lat = orgllh(1);
lon = orgllh(2);

slat = sin(lat);
clat = cos(lat);
slon = sin(lon);
clon = cos(lon);

R = [  -slat*clon  -slat*slon   clat; ...
       -slon          clon         0; ... 
       -clat*clon  -clat*slon  -slat];
        
difned = R' * ned';

ecef = difned';