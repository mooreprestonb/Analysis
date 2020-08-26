#!/usr/bin/awk -f
# write out average stdev variance skew and kurtosus from column data
BEGIN {nc = 0}
$1 != "#" {nc++;
    for(i=1;i<=NF;i++) {
	s[i] += $i;
	s2[i] += ($i)^2;
	s3[i] += ($i)^3;
	s4[i] += ($i)^4;
    }
} END { 
    print "# ",nc," Numbers: column, average, stdev, variance, skew, kurtosus";
    for (i=1;i<=NF;i++) {
	avg = s[i]/nc;
	var = s2[i]/nc - avg^2
	stdev = sqrt(var)
	skew = (s3[i]/nc - 3*avg*s2[i]/nc+2*avg^3)/(var*stdev);
	kurt = (s4[i]/nc - 4*s3[i]*avg/nc + 6*s2[i]*avg^2/nc - 3*avg^4)/(var^2)-3.
	printf "%d %f %f %f %f %f \n", i, avg, stdev, var, skew, kurt;
    }
}
