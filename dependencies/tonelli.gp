/**
 * Calculates the square root of a mod p using Tonelli's method.
 */
tonelli(a,p) = 
{
  \\ Step 0: verify that (a/p) == 1.
  if(jacobi(a,p)!=1,
		printf("tonelli: (%i/%i) == -1!", a, p);
    return(0);
  );
  
  \\ Step 1: write p-1 = 2^s * t, t odd.
  local(s = 0);
  t = p-1;
  while(t%2 == 0,
    s++;
    t = floor(t/2);
  );
  
  \\ Step 2: find b, a non-QR mod p.
  b = 2;
  while(jacobi(b,p) == 1,
    b++;
  );
    
  \\ Step 3: recurse.
  i = 2;
  c = a*b*b % p;
  for(k=1, s-1, 
    if(modexp(c, 2^(s-k-1)*t, p) == -1,
      i = i + 2^k;
      c = c*(modexp(b,2^k,p)) % p;
    );
  );
  
  \\ Step 4: calculate.
  r = (b^((i*t)/2) * a^((t+1)/2))%p;
  
  return(r);
}