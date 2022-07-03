function [x, error, iter, flag] = bicgstab(A, x, b, M, max_it, tol)

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = bicgstab(A, x, b, M, max_it, tol)
%
% bicgstab.m solves the linear system Ax=b using the 
% BiConjugate Gradient Stabilized Method with preconditioning.
%
% input   A        REAL matrix
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
%                          -1 = breakdown: rho = 0
%                          -2 = breakdown: omega = 0
%
% Updated August 2006; rbarrett@ornl.gov. (See ChangeLog for details.)
%
% =============================================================================

% ---------------
% Initialization.
% ---------------

  dim  = 0;
  iter = 0;
  flag = 0;

  alpha  = 0.0;
  beta   = 0.0;
  bnrm2  = 0.0;
  error  = 0.0;
  omega  = 1.0;
  rho    = 0.0;
  rho_1  = 0.0;

  [dim,dim] = size(A);

  p     = zeros(dim,1);
  p_hat = zeros(dim,1);
  p_tld = zeros(dim,1);
  q     = zeros(dim,1);
  r     = zeros(dim,1);
  r_tld = zeros(dim,1);
  s     = zeros(dim,1);
  s_hat = zeros(dim,1);
  t     = zeros(dim,1);
  v     = zeros(dim,1);
  z     = zeros(dim,1);
  z_tld = zeros(dim,1);

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

  r_tld = r;

  for iter = 1:max_it,

     rho   = ( r_tld'*r );
     if ( rho == 0.0 ) break, end

     % --------------------------
     % Compute direction vectors.
     % --------------------------

     if ( iter > 1 ),
        beta  = ( rho/rho_1 )*( alpha/omega );
        p = r + beta*( p - omega*v );
     else
        p = r;
     end
 
     p_hat = M \ p;
     v = A*p_hat;
     alpha = rho / ( r_tld'*v );
     s = r - alpha*v;

     % ------------------------
     % Early convergence check.
     % ------------------------

     if ( norm(s) < tol ), 
        x = x + alpha*p_hat;
        error = norm( s ) / bnrm2;
        break;
     end

     % -------------------
     % Compute stabilizer.
     % -------------------

     s_hat = M \ s;
     t = A*s_hat;
     omega = ( t'*s) / ( t'*t );

     % ---------------------
     % Update approximation.
     % ---------------------

     x = x + alpha*p_hat + omega*s_hat;

     r = s - omega*t;

     % ------------------
     % Check convergence.
     % ------------------

     error = norm( r ) / bnrm2;
     if ( error <= tol ), break, end

     if ( omega == 0.0 ), break, end
     rho_1 = rho;

  end

  if ( error <= tol | s <= tol ),

     % ----------
     % Converged.
     % ----------

     if ( s <= tol ),
        error = norm(s) / bnrm2;
     end
     flag =  0;

  elseif ( omega == 0.0 ),

     % ----------
     % Breakdown.
     % ----------

     flag = -2;

  elseif ( rho == 0.0 ),

     % ----------
     % Breakdown.
     % ----------

     flag = -1;

  else

     % ---------------
     % No convergence.
     % ---------------

     flag = 1;

  end

% --------------
% End bicgstab.m
% --------------
