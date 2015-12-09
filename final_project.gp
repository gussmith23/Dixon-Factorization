/**
 * Imports.
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
n = 135292257399511;   	


/**
 * Procedure
 */

\\ Announce!
printf("Attempting to factor %ld.\n", n);


\\\\\\\\\\\\\\\\\\\\\\\
\\ STEP 1: Choose B and M; build a factor base.

B_size = 65;
M = 5000;

\\ Expressions for Q(r) and r
Q_eqn(r) = r*r-n;
r_eqn(i) = floor(sqrt(n)) + i - M;

printf("Step 1: building a factor base. Factor base size = %i.\n", B_size);

\\ Get the factor base.
B = factor_base(n, B_size);
B_largest_prime= B[#B];
printf("Factor base of size %i:\n", B_size);
print(B);

\\\\\\\\\\\\\\\\\\\\\\\
\\ STEP 1.1: Sanity check on B and M.

\\ Via hand calculation, I calculated a such that B_largest_prime \approx n^(1/a).
\\		From this a, I looked up the probability that n factors with largest prime < n^(1/a) = B_largest_prime.
\\		This data was looked up using the Knuth/Trabb-Pardo table on p106 of Bressoud.
p_a = 3.55e-4;

\\ Now we perform the "sanity check" calculation, where we check that 
\\    expected number of factorizations >= B_size.
\\ 		This guarantees that we'll be able to Gaussian eliminate 
\\		to get zero rows.

expected_factorizations = (2*M+1) * p_a;
if(expected_factorizations < B_size, printf("Warning: sanity check on B and M failed!\n\tExpected factorizations: %d\n\tB_size:%d\n",expected_factorizations,B_size), printf("Sanity check on B and M passes!\n\tExpected factorizations: %d\n\tB_size:%d\n",expected_factorizations,B_size));


\\\\\\\\\\\\\\\\\\\\\\\\
\\ STEP 2: Quadratic Sieve.
print("Step 2: Quadratic Sieve.")
qs_sums = quadratic_sieve(B, n, M,r_eqn);
qs_threshold = (log(n)/2 + log(M)) - (3/2)*log(B_largest_prime);
printf("Quadratic Sieve produced %d sums. Attempting trial division on rows passing threshold value %f.", #qs_sums, qs_threshold);

\\ The matrix to Gaussian eliminate.
exponent_matrix = []; 

\\ The list of rs which produced viable Q(r)s.
r_list = listcreate();

\\ Iterate over each and check against threshold. If it passes, attempt trial divison.
for(i = 1, #qs_sums, {
	if(qs_sums[i] >= qs_threshold,
	
		\\ Get the number to factor.
		r = r_eqn(i);
		Qr = Q_eqn(r);
		
		\\ Create an exponent vector for it using factor base.
		local(exp_vec = vector(#B, unused, 0));
		
		\\ Check if it's negative or positive, and update exponent vector as needed.
		Qr_sign = 1;
		if(Qr < 0, 
			exp_vec[1] = 1;
			Qr_sign = -1;
		);
		
		\\ Trial division...(note that we ensure Qr is nonnegative by mult. by Qr_sign)
		local(result = trial_division( Qr_sign*Qr, B_largest_prime));
		local(td_e = result[2], td_p = result[3]);
						
		\\ If we factored completely over the factor base...
		\\ 		Note: the large primes improvement can utilize numbers that don't factor completely
		\\		over the factor base.
		if(result[4] == 1 && td_p[#td_p] <= B_largest_prime,
		
			printf("Q(%d) = %d factors completely over factor base! Placing in exponent matrix.\n", r,Qr);
			
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
			
			\\ Place in matrix.
			exponent_matrix = matconcat([exponent_matrix;exp_vec]);
			listput(r_list,r);
			
		); \\ End if factored over factor base
		
	);
});


\\ Get size.
exponent_matrix_size = matsize(exponent_matrix);

\\ IMPORTANT: first row of the matrix will be all zeros, because of how the first matconcat works. delete it.
exponent_matrix = exponent_matrix[2..exponent_matrix_size[1],1..exponent_matrix_size[2]];
exponent_matrix_size[1] -= 1; \\ make sure we update the new row count.

\\ Keep a copy of the un-reduced exp. mat.
exponent_matrix_unreduced = exponent_matrix;

printf("Quadratic Sieve produced exponent matrix with %d rows.\n", exponent_matrix_size[1]);

if(exponent_matrix_size[1] <= B_size, {
	printf("Warning: less rows than columns in the exponent matrix. When we Gaussian reduce, there is no guarantee that we'll produce zero rows.\n");
	printf("\tRows:\t%d\n", exponent_matrix_size[1]);
	printf("\tColumns:\t%d\n",B_size);
});

\\\\\\\\\\\\\\\\\\\\\\\
\\ Step 3: Gaussian Elimination.

print("Step 3: Gaussian Elimination mod 2.");

\\ Mod 2.
exponent_matrix %= 2;

\\ Append identity matrix.
exponent_matrix = matconcat([exponent_matrix, matid(exponent_matrix_size[1])]);

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
				
				exponent_matrix[j,] += exponent_matrix[found_row,];
				exponent_matrix[j,] %= 2;
								
			);
		);
		
	);	\\ End for each row.

	\\ Now we actually eliminate the row - i.e. remove it.
	if(found_row != -1,
				
			\\ CASE 1: we need to eliminate the first row.
		if(found_row == 1,
			exponent_matrix = exponent_matrix[2..exponent_matrix_size[1],];
			
			\\ CASE 2: eliminate the last row.
			, found_row == exponent_matrix_size[1],
			exponent_matrix = exponent_matrix[1..(exponent_matrix_size[1]-1),];
			
			\\ Else: eliminate a middle row.
			, exponent_matrix = matconcat([ exponent_matrix[1..(found_row-1),] ; exponent_matrix[(found_row+1)..matsize(exponent_matrix)[1],] ]);
		);
		
		if(type(exponent_matrix) == "t_POL", breakpoint());
								
		\\ Reduce row count by one.
		\\ Note that this isn't as straightforward as exponent_matrix_size = matsize(exponent_matrix)...
		\\		this is because the ACTUAL size of the exponent_matrix includes the identity portion.
		exponent_matrix_size[1] -= 1;
				
	); \\ end if foundrow != -1

}); \\ end for each column


printf("After elimination, %d rows remain.\n", exponent_matrix_size[1]);


\\\\\\\\\\\\\\\\\
\\ Step 4: Finding zero rows.

print("Step 4: Finding zero rows.");

\\ Now find zero rows.
\\ For each row.
factors = [];
for(i = 1, exponent_matrix_size[1],{
	
	\\ If it's a zero row...
	if(exponent_matrix[i,1..exponent_matrix_size[2]] == 0,
		
		\\ Generate exponent vector.
		exponents = vector(exponent_matrix_size[2], unused, 0);
		x = 1;
		y = 1;
		for(j = exponent_matrix_size[2] + 1, matsize(exponent_matrix)[2],
			if(exponent_matrix[i,j] == 1,
			
				\\ Add in the corresponding exponents
				exponents += exponent_matrix_unreduced[j-exponent_matrix_size[2],];
				
				\\ Multiply the corresponding r into x.
				x *= r_list[j-exponent_matrix_size[2]];
				x =  x % n;
			);
		);
								
		\\ Calculate value of factor base raised to exponents/2. Note that 
		\\	exponents should be divisible by 2; I'll leave the "floor" out
		\\	so that an error case where this isn't true should be immediately clear.
		exponents = exponents/2;
		
		
		for(j = 1, #B, 
			y *= B[j]^exponents[j];
			y %= n;
		);
				
		
		
		potential_factor = gcd(n,x-y);
				
		if(potential_factor != 1 && potential_factor != n,
		
			\\ If it's not already there, add it and its "complement", n/potential_factor.
			if(setsearch(factors,potential_factor) == 0, 
				factors = setunion(factors,[potential_factor]);
				factors = setunion(factors,[n/potential_factor]);			
				printf("Zero row produced factors %d and %d.", potential_factor, n/potential_factor );
			);
		);
	);
});

print("Done finding zero rows.");

\\ Sort factors.
factors = vecsort(factors);

printf("Found factors of %d:\n", n);
print(factors);
