/**
 * Knuth GCD algoritm (gcd with coefficients), useful for finding inverses.
 */
knuth_gcd(a,b) = 
{
  local(u = [1,0,a]);
  local(v = [0,1,b]);
  while(v[3] != 0,
    local(temp=v);
    v = u - floor(u[3]/v[3])*v;
    u = temp;
  );
  return(u);
}