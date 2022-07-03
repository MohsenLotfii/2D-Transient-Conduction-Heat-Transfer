function [ c, s ] = rotmat( a, b )
%
% function [ c, s ] = rotmat( a, b )
%
% rotmat.m compute the Givens rotation matrix parameters for a and b.
%
% input   a        REAL matrix element.
%         b        REAL matrix element.
%
% output  c        REAL rotation matrix element.
%         s        REAL rotation matrix element.
%
% Updated August 2006; rbarrett@ornl.gov. (See ChangeLog for details.)
% =============================================================================

   if ( b == 0.0 ),
      c = 1.0;
      s = 0.0;
   elseif ( abs(b) > abs(a) ),
      temp = a / b;
      s = 1.0 / sqrt( 1.0 + temp^2 );
      c = temp * s;
   else
      temp = b / a;
      c = 1.0 / sqrt( 1.0 + temp^2 );
      s = temp * c;
   end

% ------------
% End rotmat.m
% ------------
