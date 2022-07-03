function test_it = gmres_test()
%
% function[test_it] = gmres_test()
%
% test_qmr.m tests qmr.m. It generates several test systems and 
% applies the QMR solution algorithm as implemented in qmr.m.
%
% This and the other algorithm specific testers are modifications
% of tester.m, so some unnecessary or extraneous code is included.
%
% Created August, 2006 by Richard Barrett, rbarrett@ornl.gov.
% ============================================================================

   % ---------------
   % Initialization.
   % ---------------

   no_soln_gmres   = 0;
   guess_err_gmres = 0;                   
   num_failures       = 0;
   num_tests       = 6;

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

      if ( test == 1 | test == 4 | test == 5 | test == 6 ) 

         restrt = test*10;
         if ( restrt == 0 ) restrt = 1, end;
         [x, error, iter, flag_gmres]=gmres( A, xk, b, M, restrt, max_it, tol );
         if ( flag_gmres ~= 0 & test ~= 6 ),
            no_soln_gmres = no_soln_gmres + 1;
            'gmres failed to converge for'
            test
            num_failures = num_failures + 1;
         end
         if ( test == 6 & iter ~= 0 & flag_gmres ~= 0 )
            guess_err_gmres = guess_err_gmres + 1;
            'gmres.m  failed for initial guess = solution'
            test, iter, flag_gmres
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

   if ( no_soln_gmres ~= 0 ),
      'gmres failed test (failed to converge)',
   elseif ( guess_err_gmres ~= 0 ),
      'gmres failed test (initial guess = solution error)',
   else
      'gmres passed test';
   end

% ----------------
% End test_gmres.m
% ----------------
