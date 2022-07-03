function [x, error, iter, flag]  = jacobi(A, x, b, max_it, tol)

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
% [x, error, iter, flag]  = jacobi(A, x, b, max_it, tol)
%
% jacobi.m solves the linear system Ax=b using the Jacobi Method.
%
% input   A        REAL matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it
%
% Updated August 2006; rbarrett@ornl.gov. (See ChangeLog for details.)
%
% =============================================================================

% ---------------
% Initialization.
% ---------------

  iter = 0;
  flag = 0;
  dim = 0;

  bnrm2 = 0.0;
  error = 0.0;

  [dim,dim] = size(A);

  r   = zeros(dim,1); 
  x_1 = zeros(dim,1);
  M   = zeros(dim,1);
  N   = zeros(dim,1);

  % -----------------------------
  % Quick check of approximation.
  % -----------------------------

  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = b - A*x;
  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end

  % -----------------
  % Matrix splitting.
  % -----------------

  [ M, N ] = split( A , b, 1.0, 1 );

  % ----------------
  % Begin iteration.
  % ----------------

  for iter = 1:max_it,

     x_1 = x;

     % ---------------------
     % Update approximation.
     % ---------------------

     x   = M \ (N*x + b);

     % ------------------
     % Check convergence.
     % ------------------

     error = norm( x - x_1 ) / norm( x );
     if ( error <= tol ), break, end

  end

  % ------------------------
  % Final convergence check.
  % ------------------------

  if ( error > tol ) flag = 1; end

% ------------
% End jacobi.m
% ------------
