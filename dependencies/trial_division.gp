/**
 * Trial division on n up to a bound max_bound.
 */
trial_division(n,max_bound) =
{

	local(i=0,f=n,d=2);
	local(e=listcreate(), p=listcreate());
	
	if(n < 1, return([i,e,p,f]));
	
	if(f%d==0,
		i++;
		listput(p,d);
		listput(e,div(f,d)[1]);
		f = div(f,d)[2];		
	);
	
	d = 3;
	while(d <= max_bound && d^2 <= f,
		if(f % d == 0,
			i++;
			listput(p,d);
			listput(e,div(f,d)[1]);
			f = div(f,d)[2];		
		);
		d += 2;
	);
	
	if(f > 1 && d^2 > f,
		i++;
		listput(p,f);
		listput(e,1);
		f=1;
	);
	
	return([i,e,p,f]);

}

div(n,d) = 
{
	local(e=0,f=n);
	while(f%d == 0,
		f = f/d;
		e++;
	);
	return([e,f]);
}