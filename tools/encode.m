%  Copyright (C) 2005-2009, Axis Communications AB, LUND, SWEDEN
%
%  This file is part of Fixmath.
%
%  Fixmath is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version. You can use the comments under either
%  the terms of the GNU Lesser General Public License version 3 or
%  later, or the GNU Free Documentation License version 1.3 or later.
%
%  Fixmath is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public
%  License and a copy of the GNU Free Documentation License along
%  with Fixmath. If not, see <http://www.gnu.org/licenses/>.
%

%
% Encode polynomial coefficients in mantissa and exponent format. 
%
% encode(poly) 
%
%   poly  Polynomial coefficients obtained e.g. output from Remez' algorithm
%

function encode(vals)
     logv = log2(abs(vals));
     expn = floor(logv);
     logm = logv - expn;
     mant = 2 .^ logm;
     
     printf('\n        value       factor  sh\n')
     printf('-------------------------------\n')
 
     for k = 1:length(vals)
     
         printf('% e  ', vals(k))

         if vals(k) < 0
             printf('-')
         else
	     printf(' ')
         end
 
         m32 = round((2^31)*mant(k));
 
         printf('%#x%x  %d\n', m32 / 16, rem(m32, 16), 31 - expn(k))
     end

     printf('\n')
