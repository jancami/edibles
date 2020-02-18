def write_spectrum_ascii(output_name, x, y, yerr,header):
    f = open(output_name, 'w')
    if (header != None):
        f.write(header)
    for i in range(len(x)-1):
        if (yerr != None):
            #print x
            f.write(str('%6.4f' % x[i])+"\t"+str('%7.4f' % y[i])+"\t"+str('%7.4f' % yerr[i])+"\n")
        else:
            #print y
            f.write(str('%6.4f' % x[i])+"\t"+str('%7.4f' % y[i])+"\n")
    f.close()
    return
