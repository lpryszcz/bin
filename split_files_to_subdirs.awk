# https://stackoverflow.com/a/56062929/632242 
function first_fit(v, file) {
    # find first bin that can accomodate the volume
    for (i=1; i<=n; ++i) {
	if (b[i] > v) {
	    b[i] -= v
	    bc[i]++
	    cmd="mv "file " " i
	    print cmd
	    system(cmd)
	    return
	}
    }
    # no bin found, create new bin
    if (i > n) {
	b[++n] = c - v
	bc[n]++
	cmd="mkdir "n
	print cmd
	system(cmd)
	cmd="mv "file " " n
	print cmd
	system(cmd)
    }
        return
}
BEGIN{ if( (c+0) == 0) exit }
{ first_fit($1,$2) }
END { print "REPORT:"
    print "Created",n,"directories"
    for(i=1;i<=n;++i) print i,":", c-b[i],"bytes",bc[i],"files"
}
