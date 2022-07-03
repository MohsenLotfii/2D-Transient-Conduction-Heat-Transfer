function [x, error, iter, flag] = cgs(A, x, b, M, max_it, tol)

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = cgs(A, x, b, M, max_it, tol)
%
% cgs.m solves the linear system Ax=b using the 
% Conjugate Gradient Squared Method with preconditioning.
%
% input   A        REAL matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner
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

  A_u_hat = zeros(dim,1);
  p       = zeros(dim,1);
  p_hat   = zeros(dim,1);
  q       = zeros(dim,1);
  r       = zeros(dim,1);
  r_tld   = zeros(dim,1);
  u       = zeros(dim,1);
  u_hat   = zeros(dim,1);
  v       = zeros(dim,1);
  v_hat   = zeros(dim,1);

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

     rho = (r_tld'*r );
     if (rho == 0.0), break, end

     % --------------------------
     % Compute direction vectors.
     % --------------------------

     if ( iter > 1 ),
        beta = rho / rho_1;
        u = r + beta*q;
        p = u + beta*( q + beta*p );
     else
        u = r;
        p = u;
     end

     p_hat = M \ p;
     v_hat = A*p_hat;
     alpha = rho / ( r_tld'*v_hat );
     q = u - alpha*v_hat;
     u_hat = M \ (u+q);

     % ---------------------
     % Update approximation.
     % ---------------------

     x = x + alpha*u_hat;

     A_u_hat = A*u_hat;

     r = r - alpha*A_u_hat;

     % ------------------
     % Check convergence.
     % ------------------

     error = norm( r ) / bnrm2;
     if ( error <= tol ), break, end

     rho_1 = rho;

  end 

  if (error <= tol),

     % ----------
     % Converged.
     % ----------

     flag =  0;

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

% ---------
% End cgs.m
% ---------
