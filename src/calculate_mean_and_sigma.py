
import math, sys

#list=[1,1,0]

def main():
    sc = float(sys.argv[1])
    #list = [ float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) ]
    list = map(float, sys.argv[2:])
    #print list
    m, s = calcMeanAndSigma(list)
    if s == 0:
        print "0 variation"
    else:
        print "(%.2f - %.2f) / %.2f = %.2f" % ( sc, m, s, (sc - m) / s)

def add(x,y): return x+y

def sq(x): return x*x

def mean(list):
    n=len(list)
    return reduce(add, list)/float(n)

def sigma(list):
    n=len(list)
    m = mean(list)
    mSq=reduce(add, map(sq, list))/float(n)
    s = mSq-sq(m)
    try:
	s = math.sqrt(s)
    except:
	if s<0.0000000001:
	    s = 0
	else:
	    raise
    return s

def calc_mean_and_sigma(list):
    return mean(list), sigma(list)

if __name__ == "__main__":
    main()

