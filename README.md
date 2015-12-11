# Dixon's Factorization Algorithm with Various Improvements
MATH 467 Final Project FA15 - Prof. Brownawell

Written by Gus Smith.

Developed using Pari 2.7.4 on Windows.
Run using `gp final_project.gp`, or `\r final_project.gp` in the gp commandline.

## Description
Dixon's Factorization Algorithm including
- Pomerance's Quadratic Sieve
- Knuth/Trabb-Pardo check on factor base size
- Hensel Lifting for smaller primes

### Report
A general overview of Dixon's scheme (plus Pomerance's quadratic sieve) is as follows:
1. Create a factor base of length B of primes p such that (n/p) = 1.
2. Run the quadratic sieve over our factor base, on a total of 2M+1 values of Q(r) = r^2 - n, sqrt(n) - M < r < sqrt(n) + M.
3. Based on the results from the quadratic sieve, attempt trial division on a portion of our 2M+1 numbers. Store complete factorizations over B in an exponent matrix.
4. Reduce this matrix mod 2 to find rows of exponents representing products of Q(r)s with only even exponents.
5. Compute the square root of the product of the r-values which compose this product; call this x. Compute the square root of this even-exponent product of Q(r)s; call this y. Now gcd(n, x-y) has a chance of being a nontrivial factor of n.

In general, my program follows this procedure to a tee, with the inclusion of two main improvements:
- Knuth/Trabb-Pardo check on factor base sizes: Knuth and Trabb-Pardo produced a table of values representing probilities that a number n will factor with primes at most n^(1/a). We can use this to our advantage to check our expected number of factorizations over our factor base of size B, given 2M+1 trials. This is included in my program as a simple warning if the choices of B and M do not pass this check.
- Small primes: small primes are likely to appear multiple times in the factorization of numbers. Thus, when we run the quadratic sieve, we run it for small primes raised to higher powers. Instead of having to use the Tonelli algorithm on each higher power, we can instead use the faster method called Hensel Lifting, which will produce a "square root" mod a higher power, given a square root to the original power.


#### Report on Choice of B and M
Originally, I used Bressoud's B and M (30, 5000). I didn't expect this to work for my number; however, to my surprise, it does. After Gaussian elimination, there is a single zero row, which luckily produces one of our two factors. The "sanity check" on B and M fails miserably, of course, predicting less than one complete factorization over the factor base. However, because this choice of B and M works, I decided to stick with it.

Larger choices of B and M work equally well, but take much more time, and have the potential to use a lot of memory in the Gaussian reduction phase.

#### Potential Improvements
There are many potential improvements to my own implementation. Each subroutine has room for optimization; especially the quadratic sieve. I was more focused on producing a working factorization than a working AND optimal factorization - premature optimization is a dangerous thing!

A serious improvement which I nearly implemented was in the Gaussian elimination phase. Instead of appending an identity matrix to the exponent matrix, I considered using lists to indicate which rows had contributed their sum to a given row. For a large exponent matrix, an appended identity matrix is very sparse, and can take up a lot of unused room. Lists as I described performed up to an order of magnitude better on space, from some initial testing I did. 

### Conclusion
It was remarkable to see all of our knowledge come together like this in the end. The first time I had an accurate factorization, I felt overjoyed. It is always satisfying when the things learned in a course can be combined into a tangible result!
