/**
 * Imports.
 * These may need to be done in a specific order.
 */
\r dependencies/modexp.gp;
\r dependencies/prime_sieve.gp;
\r dependencies/hensel_lifting.gp;
\r dependencies/jacobi.gp;
\r dependencies/knuth_gcd.gp;
\r dependencies/factor_base.gp;
\r dependencies/tonelli.gp;
\r dependencies/quadratic_sieve.gp;
\r dependencies/trial_division.gp;
\r dependencies/helpers.gp;


/**
 * Fields
 */

\\ The number to factor.
\\n = 135292257399511;   	\\ Assigned number.
n = 499 94860 12441;			\\ Bressoud's number (Fact. and Primality Testing p110)
\\n = 100;



/**
 * Procedure
 */

\\ Announce!
printf("Attempting to factor %ld.", n);


\\\\\\\\\\\\\\\\\\\\\\\
\\ STEP 1: Choose B and M; build a factor base.
\\B_size = 10; 		\\ TODO figure out how to work the trabb-pardo/knuth table into this.
\\M = 20;						\\ These are bad choices i'm sure; just dummies for now.

\\ Bressoud's B and M (p110)
B_size = 30;
M = 5000;

\\ Expressions for Q(r) and r
Q_eqn(r) = r*r-n;
r_eqn(i) = floor(sqrt(n)) + i - M;

printf("Step 1: building a factor base. Factor base size = %i.", B_size);

\\ Get the factor base.
B = factor_base(n, B_size);
B_largest_prime= B[#B];
printf("Factor base of size %i:", B_size);
print(B);

\\\\\\\\\\\\\\\\\\\\\\\
\\ STEP 1.1: Sanity check on B and M.


\\\\\\\\\\\\\\\\\\\\\\\\
\\ STEP 2: Quadratic Sieve.
qs_sums = quadratic_sieve(B, n, M,r_eqn);
qs_threshold = (log(n)/2 + log(M)) - (3/2)*log(B_largest_prime);


\\ The matrix to Gaussian eliminate.
exponent_matrix = []; 
Qr_list = []; \\ The list of Q(r)s that go with each exponent mat.
r_list = [];
i_list = [];

factored_completely = 0; \\ Debug: count the number of Q(r)s that factor completely.

\\ Iterate over each and check against threshold. If it passes, attempt trial divison.
for(i = 1, #qs_sums, {
	if(qs_sums[i] >= qs_threshold,
	
		\\ Get the number to factor.
		r = r_eqn(i);
		Qr = Q_eqn(r);
		printf("Attempting trial division on Q(%d) = %d\n", r,Qr);
		
		
		\\ Create an exponent vector for it using factor base.
		local(exp_vec = vector(#B, unused, 0));
		
		\\ Check if it's negative or positive, and update factor base as needed.
		Qr_sign = 1;
		if(Qr < 0, 
			exp_vec[1] = 1;
			Qr_sign = -1;
		);
		
		\\ Trial division...(note that we ensure Qr is nonnegative by mult. by Qr_sign)
		\\printf("Attempting trial division on %i...\n", Qr_sign*Qr);
		local(result = trial_division( Qr_sign*Qr, B_largest_prime));
		local(td_e = result[2], td_p = result[3]);
		
		print(result);
				
		\\ If we factored completely over the factor base.
		\\ TODO is this right?
		if(result[4] == 1 && td_p[#td_p] <= B_largest_prime,
		
			factored_completely += 1;
		
			\\ Now we put the exponents from trial division into our exponent matrix.
			local(i_tde = 1, i_expvec = 1);
			while(i_tde <= #td_e && i_expvec <= #exp_vec,
				if(td_p[i_tde] == B[i_expvec],
					\\ Case 1: The prime in the trial division equals the prime in our factor base.
					exp_vec[i_expvec] = td_e[i_tde];
					i_expvec++; i_tde++;,
					\\ ELSE If
					\\ Case 2
					td_p[i_tde] < B[i_expvec],
					i_tde++,
					\\ ELSE If
					\\ Case 3
					td_p[i_tde] > B[i_expvec],
					i_expvec++
				);
			);
			
			\\print("Exponent vector to be placed in table:");
			\\print(exp_vec);	
			
			\\ Place in matrix.
			exponent_matrix = matconcat([exponent_matrix;exp_vec]);
			Qr_list = matconcat([Qr_list,Qr]);
			r_list = matconcat([r_list,r]);
			i_list = matconcat([i_list,i]);
			
			
		); \\ End if factored over factor base
		
	);
});


\\ Get size.
exponent_matrix_size = matsize(exponent_matrix);

\\ IMPORTANT: first row of the matrix will be all zeros, because of how the first matconcat works. delete it.
exponent_matrix = exponent_matrix[2..exponent_matrix_size[1],1..exponent_matrix_size[2]];
exponent_matrix_size[1] -= 1; \\ make sure we update the new row count.

\\ Additionally, change qr to a vector, not a 1xk matrix.
Qr_list = Qr_list[1,];
r_list = r_list[1,];
i_list = i_list[1,];

\\ Keep a copy of the un-reduced exp. mat.
exponent_matrix_unreduced = exponent_matrix;

\\ Mod 2.
exponent_matrix %= 2;


\\ OBSOLETE:
\\ Append identity matrix.
\\exponent_matrix = matconcat([exponent_matrix, matid(exponent_matrix_size[1])]);

\\ New strategy: Keep a vector of lists which will do the same as the above identity
\\		matrix; that is, track what rows are added together to produce this row.
\\ 		This will (hopefully) reduce memory consumption.
component_lists = vector(exponent_matrix_size[1], unused, listcreate());
for(i = 1, exponent_matrix_size[1], listput(component_lists[i],i));

\\ Eliminate mod 2.
\\ For each column (starting with larger primes first)...
forstep(i = #B, 1, -1,{
	
	found_row = -1;
	
	\\ For each row...
	for(j = 1, exponent_matrix_size[1],
	
		\\ If we need to eliminate...
		if(exponent_matrix[j, i] == 1,
		
			\\ Check whether or not an above row has a "1" in this column.
			\\		If there is such a row, we'll use that row to eliminate this one.
			\\		Otherwise, we'll log this row in found_row and use it to eliminate later rows.
			if(found_row == -1, 

				found_row = j;
				
				, \\ <-- begin ELSE clause. (i don't like this syntax!)
				\\ ELSE: found_row equals some j, and we should use that row j to eliminate this row.
				
				\\ For each entry in row j, add the corresponding entry from found_row, mod 2.
				for(k = 1, matsize(exponent_matrix)[2], 
					exponent_matrix[j,k] += exponent_matrix[found_row,k];
					exponent_matrix[j,k] %= 2;
				);
				
				\\ Log that we added row found_row into row j.
				listput(component_lists[j], found_row);
				
			);
		);
	);
});


\\ Now find zero rows.
\\ For each row.
factors = [];
for(i = 1, exponent_matrix_size[1],{
	
	\\ If it's a zero row...
	if(exponent_matrix[i,1..exponent_matrix_size[2]] == 0,
	
		print("zero row");
		
		\\ Generate exponent vector.
		exponents = vector(exponent_matrix_size[2], unused, 0);
		
		\\ The list of rows that went into making this row.
		component_list = component_lists[i];
		
		\\ For each entry in component_list...
		for(j = 1, #component_list, 
			exponents += exponent_matrix_unreduced[component_list[j],];
		);
						
		\\ Calculate value of factor base raised to exponents/2. Note that 
		\\	exponents should be divisible by 2; I'll leave the "floor" out
		\\	so that an error case where this isn't true should be immediately clear.
		exponents = exponents/2;
		
		x = 1;
		for(j = 1, #B, 
			x *= B[j]^exponents[j];
		);
		
		y = r_list[i];
		
		potential_factor = gcd(n,x-y);
		
		if(potential_factor != 1 && potential_factor != n,
		
			if(setsearch(factors,potential_factor) == 0, factors = setunion(factors,[potential_factor]));
			
		);
	);

});


\\ Sort factors.
factors = vecsort(factors);

print("Found factors:");
print(factors);
