function test_it = qmr_test()
%
% function[test_it] = qmr_test()
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

   no_soln_qmr   = 0;
   guess_err_qmr = 0;                   
   num_failures     = 0;
   num_tests     = 6;

   ep  = eps;
   tol = ep * 1000;

   % -------------------
   % Form linear system.
   % -------------------

   for test = 1:num_tests

      if ( test == 1 | test == 4 | test == 5 | test == 6 ),

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
   
         % ---------
         % Test QMR.
         % ---------
   
         [x, error, iter, flag_qmr] = qmr(A, xk, b, M, max_it, tol);
         if ( flag_qmr ~= 0 & test ~= 6 ),
            no_soln_qmr = no_soln_qmr + 1;
            'qmr failed to converge for'
            test, error
            num_failures = num_failures + 1;
         end
         if ( test == 6 & iter ~= 0 & flag_qmr ~= 0 ),
            guess_err_qmr = guess_err_qmr + 1;
            'qmr.m failed for initial guess = solution'
            test, iter, flag_qmr
            num_failures = num_failures + 1;
         end
      end
   end
 
   % ------------------------
   % Write outcome to screen.
   % ------------------------

   if ( no_soln_qmr ~= 0 ),
      'QMR failed test (failed to converge)',
   elseif ( guess_err_qmr ~= 0 ),
      'qmr.m failed test (initial guess = solution error)',
   else
      'QMR passed test'
   end

% -------------
% End qmr_test.
% -------------
