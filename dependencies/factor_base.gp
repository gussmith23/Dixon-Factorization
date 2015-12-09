/**
 * Builds a factor base of B primes p such that (n/p) = 1.
 */
factor_base(n,B) =
{
	\\ The max prime for the prime sieve to find.
	lim = 0;

	\\ Our list of primes returned by the prime sieve.
	local(h = []);
	local(last_prime_checked = 0);

	\\ The total number of primes p such that (n/p) == 1 which we are searching for.
	\\ Note that we also add -1 into our factor base.
	local(goal = B);

	\\ This list will contain the primes p such that (n/p) == 1.
	local(factor_base = listcreate(goal));
	\\ Improvement: include -1 in factor base.
	listput(factor_base, -1);
	\\ TODO: should the factor base include 2, too?

	while (#factor_base<goal,

		\\ Find more primes
		lim += 100;
		h = prime_sieve(lim);
		
		\\ Check (n/p) on new primes
		for(i = last_prime_checked+1, #h,
			if(jacobi(n,h[i])==1, 
				listput(factor_base,h[i]);
				if(#factor_base >= goal, break());
			);
		);
		
		\\ This value makes sure we don't check any primes we've already checked
		last_prime_checked = #h;

	);

	\\ Output.
	\\print("Factor base:");
	\\print(factor_base);
	
	return(factor_base);


	/**
	 * Equivalent code using Pari functions:
	 *
	test = listcreate(100);
	forprime(p=2,100000,{
		if(kronecker(n,p) == 1,
			listput(test,p);
			if(#test == 100, break());
		);
	});
	print(test); */
}