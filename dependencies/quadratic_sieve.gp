hensel_lift_amounts = vector(23,unused_var,1);
hensel_lift_amounts[2] = 7;
hensel_lift_amounts[3] = 4;
hensel_lift_amounts[5] = 3;
hensel_lift_amounts[7] = 2;
hensel_lift_amounts[11] = 2;
hensel_lift_amounts[13] = 2;
hensel_lift_amounts[17] = 2;
hensel_lift_amounts[19] = 2;
hensel_lift_amounts[23] = 2;

/**
 * Quadratic sieve over factor_base for r(x) in the range sqrt(n) - M <= r(x) <= sqrt(n) + M.
 */
 quadratic_sieve(factor_base,n,M,r) =
{
	local(qs_mat = matrix(2*M+1, #factor_base));
	local(sqrt_n = floor(sqrt(n)));
	
	for(i=2, #factor_base, 		\\ Note: we skip -1 by starting at index 2.
	
		\\ The prime for this iteration.
		local(p = factor_base[i]);
		
		\\ The number of times the prime should be hensel lifted.
		local(hensel_lift_amount); \\ This number actually specifies how many times the following loop should run.
		if (p > 23, hensel_lift_amount = 1, hensel_lift_amount = hensel_lift_amounts[p]);
		
		
		\\ Our solution.
		local(a);
		
		for(j = 1, hensel_lift_amount,		\\ Note that this loop runs only once for primes greater than 23.
			
			\\ We want to find an a such that n = a^2 mod p. We know this a exists due to how we 
			\\		picked our factor base.
			\\ There are two methods for calculating a:
			\\ METHOD 1: Calculate directly using Tonelli's method.
			\\ METHOD 2: if we have an a mod p^k, we can find a solution mod p^{k+1} via Hensel lifting.
						
			if (j==1,
				\\ METHOD 1. We use this on the first iteration.
				a = tonelli(n,p);
				,
				\\ else METHOD 2. 
				a = hensel_lift(a,n,p,j-1);
			);
			
			\\ Now, find all rows whose corresponding r is equivalent to a.
			\\ r = sqrt_n + m - M - 1 
			\\ TODO don't do it this way... find the first r and then do forstep.
			
			val_to_add = ceil(log(p));
			
			for(m = 1,2*M+1,
				\\ Check for equivalence with a mod p.
				if(((r(m))%p) == a, qs_mat[m,i] += val_to_add);
				\\ Check for equivalence with -a mod p.
				if(((r(m))%p) == p-a, qs_mat[m,i] += val_to_add);
			);
		);
		
	);
	
	return(sum(x = 1, matsize(qs_mat)[2], qs_mat[1..2*M+1,x]));
}