function [x, error, iter, flag] = qmr( A, x, b, M, max_it, tol )

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = qmr( A, x, b, M, max_it, tol )
%
% qmr.m solves the linear system Ax=b using the 
% Quasi Minimal Residual Method with preconditioning.
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
%         flag     INTEGER: 0: solution found to tolerance
%                           1: no convergence given max_it
%                     breakdown:
%                          -1: rho
%                          -2: beta
%                          -3: gamma
%                          -4: delta
%                          -5: ep
%                          -6: xi
%
% Updated August 2006; rbarrett@ornl.gov. (See ChangeLog for details.)
%
% =============================================================================

% ---------------
% Initialization.
% ---------------

   dim  = 0;
   flag = 0;
   iter = 0;

   beta         = 0.0;
   bnrm2        = 0.0;
   delta        = 0.0;
   ep           = 0.0;
   error        = 0.0;
   eta          = 0.0;
   gamma        = 0.0;
   gamma_1      = 0.0;
   rho          = 0.0;
   rho_1        = 0.0;
   theta        = 0.0;
   theta_1      = 0.0;
   xi           = 0.0;

   [dim,dim] = size(A);
 
   d            = zeros(dim,1);
   p            = zeros(dim,1);
   p_tld        = zeros(dim,1);
   q            = zeros(dim,1);
   r            = zeros(dim,1);
   s            = zeros(dim,1);
   v            = zeros(dim,1);
   v_tld        = zeros(dim,1);
   w            = zeros(dim,1);
   w_tld        = zeros(dim,1);
   y            = zeros(dim,1);
   y_tld        = zeros(dim,1);
   z            = zeros(dim,1);
   z_tld        = zeros(dim,1);
 
   % -----------------------------
   % Quick check of approximation.
   % -----------------------------

   bnrm2 = norm( b );
   if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

   r = b - A*x;
   error = norm( r ) / bnrm2;
   if ( error < tol ) return, end

   % ---------------------------------------------------
   % Factor preconditioner into left and right matrices.
   % ---------------------------------------------------

   [M1,M2] = lu( M );

   v_tld = r;
   y     = M1 \ v_tld;
   rho   = norm( y );

   w_tld = r;
   z     = M2' \ w_tld;
   xi    = norm( z );

   gamma =  1.0;
   eta   = -1.0;
   theta =  0.0;

  % ----------------
  % Begin iteration.
  % ----------------

   for iter = 1:max_it,

      if ( rho == 0.0 | xi == 0.0 ), break, end

      v = v_tld / rho;
      y = y / rho;

      w = w_tld / xi;
      z = z / xi;

      delta = z'*y;
      if ( delta == 0.0 ), break, end

      y_tld = M2 \ y;
      z_tld = M1'\ z;

     % --------------------------
     % Compute direction vectors.
     % --------------------------

      if ( iter > 1 ),
         p = y_tld - ( xi*delta / ep )*p;
         q = z_tld - ( rho*delta / ep )*q;
      else
         p = y_tld;
         q = z_tld;
      end

     % --------------------
     % Check for breakdown.
     % --------------------

      p_tld = A*p;
      ep = q'*p_tld;
      if ( ep == 0.0 ), break, end

      beta = ep / delta;
      if ( beta == 0.0 ), break, end

      v_tld = p_tld - beta*v;
      y =  M1 \ v_tld;

      rho_1 = rho;
      rho = norm( y );
      w_tld = ( A'*q ) - ( beta*w );
      z =  M2' \ w_tld;

      xi = norm( z );

      gamma_1 = gamma;
      theta_1 = theta;

      theta = rho / ( gamma_1*beta );
      gamma = 1.0 / sqrt( 1.0 + (theta^2) );
      if ( gamma == 0.0 ), break, end

      eta = -eta*rho_1*(gamma^2) / ( beta*(gamma_1^2) );

      % ---------------------------
      % Compute adjustment vectors.
      % ---------------------------

      if ( iter > 1 ),
         d = eta*p + (( theta_1*gamma )^2)*d;
         s = eta*p_tld + (( theta_1*gamma )^2)*s;
      else
         d = eta*p;
         s = eta*p_tld;
      end

     % ---------------------
     % Update approximation.
     % ---------------------

      x = x + d;

     % ------------------
     % Check convergence.
     % ------------------

      r = r - s;
      error = norm( r ) / bnrm2;
      if ( error <= tol ), break, end

   end

   % ---------------------------
   % Set iteration outcome flag. 
   % ---------------------------

   if ( error <= tol ),                        % converged
      flag =  0;
   elseif ( rho == 0.0 ),                      % breakdown
      flag = -1;
   elseif ( beta == 0.0 ),
      flag = -2;
   elseif ( gamma == 0.0 ),
      flag = -3;
   elseif ( delta == 0.0 ),
      flag = -4;
   elseif ( ep == 0.0 ),
      flag = -5;
   elseif ( xi == 0.0 ),
      flag = -6;
   else                                        % no convergence
      flag = 1;
   end

% ----------
%  End qmr.m
% ----------
