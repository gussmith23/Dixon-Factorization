# Dixon's Factorization Algorithm with Various Improvements
MATH 467 Final Project FA15 - Prof. Brownawell

Developed using Pari 2.7.4 on Windows.
Run using `gp final_project.gp`, or `\r final_project.gp` in the gp commandline.

## Description
Dixon's Factorization Algorithm including
- Pomerance's Quadratic Sieve
- Knuth/Trabb-Pardo check on factor base size

### Choice of B and M
Originally, I used Bressoud's B and M (30, 5000). I didn't expect this to work for my number; however, to my surprise, it does. After Gaussian elimination, there is a single zero row, which luckily produces one of our two factors. The "sanity check" on B and M fails miserably, of course, predicting less than one complete factorization over the factor base. However, because this choice of B and M works, I decided to stick with it.
