/**
 * Fast modular exponentiaton.
 */
modexp(b,e,m) =
{
  local(n = 1);
  while(e != 0,
    if(e%2 == 1,
      n = (n*b)%m;
    );
    e = floor(e/2);
    b = (b*b)%m;
  );
  return(n);
}