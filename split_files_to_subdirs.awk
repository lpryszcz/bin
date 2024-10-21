# https://stackoverflow.com/a/56062929/632242
header = "#!\n[operation]\naction=archive\ncompression=no\n[objects]"

function first_fit(v, file) {
    # find first bin that can accomodate the volume
    for (i=1; i<=n; ++i) {
	if (b[i] > v) {
	    b[i] -= v
	    bc[i]++
	    cmd = "echo "file" >> _toTape."i".sh"
	    system(cmd)#; print cmd
	    return
	}
    }
    # no bin found, create new bin
    if (i > n) {
	b[++n] = c - v
	bc[n]++
	cmd = "echo '"header"' > _toTape."i".sh"
	system(cmd)#; print cmd
	cmd = "echo "file" >> _toTape."i".sh"
	system(cmd)#; print cmd
    }
        return
}
BEGIN{ if( (c+0) == 0) exit }
{ first_fit($1,$2) }
END { print "REPORT:"
    print "Created", n, "jobs"
    for(i=1;i<=n;++i) print i, ":", c-b[i], "bytes", bc[i], "path(s)"
}
