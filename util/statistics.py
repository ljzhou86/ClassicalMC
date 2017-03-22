import numpy

def statistics(index):
    E =[] # energy
    E2=[] # energy^2
    M =[] # magnetization
    M2=[] # magnetization^2
    beta=0.0
    temperature=0.0
    #print index
    for line in open("output.%d" % (index) ):

        if "temperature" in line:
            tmp = line.strip().split()
            beta = float( tmp[-2] )
            temperature = float( tmp[3] )

        tmp = line.strip().split()
        try:
            e = float( tmp[0] )
            m = float( tmp[1] ) * float( tmp[1] )
            m+= float( tmp[2] ) * float( tmp[2] )
            m+= float( tmp[3] ) * float( tmp[3] )
            E.append(e)
            E2.append(e*e)
            M.append(m)
            M2.append(m*m)
        except:
            pass

    M_avg = numpy.mean(M)
    Cv = beta*beta*( numpy.mean(E2) - numpy.mean(E)*numpy.mean(E) )
    Chi= beta*( numpy.mean(M2) - numpy.mean(M)*numpy.mean(M) )
    return temperature, M_avg, Cv, Chi

if __name__ == '__main__':
    import sys
    for i in sys.argv[1:]:
        temperature, M_avg, Cv,Chi = statistics( int(i) )
        print i,temperature, M_avg, Cv, Chi
