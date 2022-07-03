function [ M, N, b ] = split( A, b, w, flag )
%
% function [ M, N, b ] = split( A, b, w, flag )
%
% split.m sets up the matrix splitting for the stationary
% iterative methods: Jacobi and SOR (Gauss-Seidel when w = 1.0 )
%
% input   A        REAL matrix
%         b        REAL right hand side vector (for SOR)
%         w        REAL relaxation scalar
%         flag     INTEGER flag for method: 1 = jacobi
%                                           2 = sor
%
% output  M        REAL matrix
%         N        REAL matrix such that A = M - N
%         b        REAL rhs vector ( altered for SOR )
%
% Updated August 2006; rbarrett@ornl.gov. (See ChangeLog for details.)
% =============================================================================

  if ( flag == 1 ),

     % -----------------
     % Jacobi splitting.
     % -----------------

     M = diag(diag(A));
     N = diag(diag(A)) - A;

  elseif ( flag == 2 ),

     % ---------------------------
     % SOR/Gauss-Seidel splitting.
     % ---------------------------

     b = w * b;
     M =  w * tril( A, -1 ) + diag(diag( A ));
     N = -w * triu( A,  1 ) + ( 1.0 - w ) * diag(diag( A ));

  end;

% -----------
% End split.m
% -----------
