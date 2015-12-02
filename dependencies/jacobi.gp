/**
 * Returns the Jacobi symbol (a/m).
 */ 
jacobi(a,m) =
{
  local(a = a%m);
  local(t = 1);
  while(a!=0,
    local(c = 0);
    while(a%2 == 0,
      a = a/2;
      c = 1-c;
    );
    if(c==1,
      if( (m%8 == 3) || (m%8 == 5) , t = -t)
    );
    if( (a%4 == m%4) && (a%4 == 3) , t = -t);
    
    local(temp = m);
    m = a;
    a = temp;
    a = a%m;
  );
  
  local(r = 0);
  if(m==1,r = t);
  return(r);
}