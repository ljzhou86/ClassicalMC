def autocorrelation(O,j):
    O_avg = 0.0
    O2_avg = 0.0
    OO_avg = 0.0
    n=len(O)
    for i in xrange(n):
        O_avg  += O[i]/n
        O2_avg += O[i]*O[i]/n
    for i in xrange(n-j-1):
        OO_avg += O[i]*O[i+j]/(n-j-1)

    return (OO_avg - O_avg * O_avg) / ( O2_avg - O_avg * O_avg)

if __name__ == '__main__':
    import sys
    import math
    f=open(sys.argv[1],'r')

    E=[]
    M=[]
    for line in f:
        try:
            tmp=line.strip().split()
            E.append( float(tmp[0]) )
            m = float(tmp[1])*float(tmp[1]) + float(tmp[2])*float(tmp[2]) + float(tmp[3])*float(tmp[3])
            M.append( math.sqrt(m) )
        except:
            pass

        
    for j in xrange(100):
        corrE = autocorrelation(E,j)
        corrM = autocorrelation(M,j)
        print j, corrE, corrM
