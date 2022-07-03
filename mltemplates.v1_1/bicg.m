function [x, error, iter, flag] = bicg(A, x, b, M, max_it, tol)

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = bicg(A, x, b, M, max_it, tol)
%
% bicg.m solves the linear system Ax=b using the 
% BiConjugate Gradient Method with preconditioning.
%
% input   A        REAL matrix
%         M        REAL preconditioner matrix
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
%                          -1 = breakdown
%
% Updated August 2006; rbarrett@ornl.gov. (See ChangeLog for details.)
% =============================================================================

% ----------------
%  Initialization.
% ----------------

  dim  = 0;
  iter = 0;
  flag = 0;

  alpha = 0.0;
  beta  = 0.0;
  bnrm2 = 0.0;
  error = 0.0;
  rho   = 0.0;
  rho_1 = 0.0;

  [dim,dim] = size(A);

  r     = zeros(dim,1);
  r_tld = zeros(dim,1);
  p     = zeros(dim,1);
  p_tld = zeros(dim,1);
  q     = zeros(dim,1);
  z     = zeros(dim,1);
  z_tld = zeros(dim,1);

  % -----------------------------
  % Quick check of approximation.
  % -----------------------------

  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = b - A*x;

  error = norm( r ) / bnrm2;
  if ( error < tol ), return, end

  r_tld = r;

  % ----------------
  % Begin iteration.
  % ----------------

  for iter = 1:max_it

     z = M \ r;
     z_tld = M' \ r_tld;
     rho   = ( z'*r_tld );
     if ( rho == 0.0 ),
        iter, error, 'rho breakdown'
        break
     end

     % --------------------------
     % Compute direction vectors.
     % --------------------------

     if ( iter > 1 ),
        beta = rho / rho_1;
        p   = z  + beta*p;
        p_tld = z_tld + beta*p_tld;
     else
        p = z;
        p_tld = z_tld;
     end

     % -----------------------
     % Compute  residual pair.
     % -----------------------

     q = A*p;
     q_tld = A'*p_tld;
     alpha = rho / (p_tld'*q );

     % ---------------------
     % Update approximation.
     % ---------------------

     x = x + alpha*p;
     r = r - alpha*q;
     r_tld = r_tld - alpha*q_tld;

     % ------------------
     % Check convergence.
     % ------------------

     error = norm( r ) / bnrm2;
     if ( error <= tol ), break, end

     rho_1 = rho;

  end

  if ( error <= tol ),

     % ----------
     % Converged.
     % ----------

     flag =  0;

  elseif ( rho == 0.0 ),

     % ----------
     % Breakdown.
     % ----------

     flag = -1;

  else,

     % ---------------
     % No convergence.
     % ---------------

     flag = 1;

  end

% ----------
% End bicg.m
% ----------
