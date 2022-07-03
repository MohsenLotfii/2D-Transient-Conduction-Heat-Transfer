function test_it = test_jacobi()
%
% function[test_it] = test_jacobi()
%
% test_jacobi.m tests jacobi.m. It generates several test systems and 
% applies the Jacbi solution algorithm as implemented in jacobi.m.
%
% This and the other algorithm specific testers are modifications
% of tester.m, so some unnecessary or extraneous code is included.
%
% Created August, 2006 by Richard Barrett, rbarrett@ornl.gov.
% ============================================================================

   % ---------------
   % Initialization.
   % ---------------

   no_soln_jac   = 0;                   % convergence check
   guess_err_jac = 0;                   % initial guess = solution error
   num_failures     = 0;
   num_tests     = 6;

   ep  = eps;
   tol = ep * 1000;

   % -----------------------------------------------------------------
   % Iterate over linear systems, applying solution algorithm to each.
   % -----------------------------------------------------------------

   for test = 1:num_tests
      test
      A = matgen( test*10 );             % form test matrix
      [sizeA,sizeA] = size(A);
      max_it = sizeA * 10;
      normA = norm( A,inf );
      if ( test == 1 | test == 2 | test == 3 | test == 6 ),
         for i = 1:sizeA,                % set rhs = row sums
            temp = 0.0;
            for j = 1:sizeA,
               temp = temp + A(i,j);
            end
            b(i,1) = temp;
         end
      else 
         b = ones(sizeA,1);              % set rhs = unit vector
      end
      if ( test < 4 ),
         M = eye(sizeA);                 % no preconditioning
      else
         M = diag(diag(A));              % diagonal preconditioning
      end

      if ( test < 6 ),
         xk = zeros(sizeA,1);            % initial guess = zero vector
      else
         xk = A \ b;                     % initial guess = solution
      end

      % ----------------
      % Apply algorithm.
      % ----------------

      if ( test == 1 | test == 6 ),
         jac_maxit = max_it * 3;
         [x, error, iter, flag_jac] = jacobi(A, xk, b, jac_maxit, tol);
         if ( flag_jac ~= 0 & test ~= 6 ),
            no_soln_jac = no_soln_jac + 1;
            'jacobi failed to converge for'
            test, error
            num_failures = num_failures + 1;
         end
         if ( test == 6 & iter ~= 0 & flag_jac ~= 0 )
            guess_err_jac = guess_err_jac + 1;
            'jacobi.m  failed for initial guess = solution'
            test, iter, flag_jac
            num_failures = num_failures + 1;
         end
      end
   end

   % -----------------------
   % Print results to stdio.
   % -----------------------

   TESTING = '             COMPLETE'

   if ( num_failures == 0 ),
      RESULTS = '             ALL TESTS PASSED', end

   if ( no_soln_jac ~= 0 ),
      'jacobi failed test (failed to converge)',
   elseif ( guess_err_jac ~= 0 ),
      'jacobi failed test (initial guess = solution error)',
   else
      'jacobi passed test';
   end

% -----------------
% End test_jacobi.m
% -----------------
