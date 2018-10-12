## Copyright (C) 2018 juho
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Author: juho <iipponen@princeton.edu>
## Created: 2018-08-10

function [streamfunction] = solveZonalStreamfunction (u_div, omega, lat_vec, lon_vec, p_vec, lat1, lat2, walkerFromZonalPert)


% Subtract zonal means
if walkerFromZonalPert
  u_div = u_div - mean(u_div,1);
  omega = omega - mean(omega,1);
end




end
