/**
 * Hensel lifting.
 * input: t, a, p, k such that
 * t^2 \equiv a \mod p^k.
 * returns x such that x^2 \equiv a \mod p^{k+1}.
 */
hensel_lift(t,a,p,k) = 
{
  local(l);
  local(inverse);
  inverse = (knuth_gcd(2*t, p)[1]) % p;
  l = (inverse*(a-t^2)/(p^k));
  return( (t+l*p^k)%p^(k+1) );
}