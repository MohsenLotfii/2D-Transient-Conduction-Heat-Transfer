function test_it = test_cgs()
%
% function[test_it] = test_cgs()
%
% test_cgs.m tests cgs.m. It generates several test systems and 
% applies the CGS solution algorithm as implemented in cgs.m.
%
% This and the other algorithm specific testers are modifications
% of tester.m, so some unnecessary or extraneous code is included.
%
% Created August, 2006 by Richard Barrett, rbarrett@ornl.gov.
% ============================================================================

   % ---------------
   % Initialization.
   % ---------------

   no_soln_cgs   = 0;
   guess_err_cgs = 0;                   
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

      if ( test == 1 | test == 4 | test == 5 | test == 6 ) 

         [x, error, iter, flag_cgs] = cgs(A, xk, b, M, max_it, tol);
         if ( flag_cgs ~= 0 & test ~= 6 ),
            no_soln_cgs = no_soln_cgs + 1;
            'cgs failed to converge for'
            test, error
            num_failures = num_failures + 1;
         end
         if ( test == 6 & iter ~= 0 & flag_cgs ~= 0 )
            guess_err_cgs = guess_err_cgs + 1;
            'cgs.m  failed for initial guess = solution'
            test, iter, flag_cgs
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

   if ( no_soln_cgs ~= 0 ),
      'cgs failed test (failed to converge)',
   elseif ( guess_err_cgs ~= 0 ),
      'cgs failed test (initial guess = solution error)',
   else
      'cgs passed test';
   end

% --------------
% End test_cgs.m
% --------------
