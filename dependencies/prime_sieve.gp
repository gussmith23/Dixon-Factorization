/**
 * Returns a list of primes up to bound n.
 */
prime_sieve(n) =
{
  v=vector(n, i, 1);
  v[1]=0;
  forprime(p=2, sqrt(n),
    forstep(j=p+p,n,p,v[j]=0)
  );
  h=listcreate(n);
  for(i=1,n,if(v[i],listput(h,i)));
  return(h);
}