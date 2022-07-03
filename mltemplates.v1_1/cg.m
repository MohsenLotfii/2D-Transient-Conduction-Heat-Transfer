function [x, error, iter, flag] = cg(A, x, b, M, max_it, tol)

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = cg(A, x, b, M, max_it, tol)
%
% cg.m solves the symmetric positive definite linear system Ax=b 
% using the Conjugate Gradient method with preconditioning.
%
% input   A        REAL symmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
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
% Initializations
% ---------------

  dim = 0;
  flag = 0;
  iter = 0;

  alpha = 0.0;
  beta  = 0.0;
  bnrm2 = 0.0;
  error = 0.0;
  rho   = 0.0;
  rho_1 = 0.0;

  [dim,dim] = size(A);

  p     = zeros(dim,1);
  q     = zeros(dim,1);
  r     = zeros(dim,1);
  z     = zeros(dim,1);

  % -----------------------------
  % Quick check of approximation.
  % -----------------------------

  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = b - A*x;
  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end

  % ----------------
  % Begin iteration.
  % ----------------

  for iter = 1:max_it

     z  = M \ r;
     rho = (r'*z);

     % -------------------------
     % Compute direction vector.
     % -------------------------

     if ( iter > 1 ),
        beta = rho / rho_1;
        p = z + beta*p;
     else
        p = z;
     end

     q = A*p;
     alpha = rho / (p'*q );

     % ---------------------
     % Update approximation.
     % ---------------------

     x = x + alpha * p;

     r = r - alpha*q;

     % ------------------
     % Check convergence.
     % ------------------

     error = norm( r ) / bnrm2;
     if ( error <= tol ), break, end 

     rho_1 = rho;

  end

  % ------------------------
  % Final convergence check.
  % ------------------------

  if ( error > tol ) flag = 1; end

% --------
% End cg.m
% --------
