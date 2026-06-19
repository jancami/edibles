from importdata import importdata
def try_orders(target,minrange,maxrange,ContinuumMin,ContinuumMax):
    for x in range(24):
        try:
            importdata(target=target, minrange=minrange, maxrange=maxrange,
                                                                                  data_piece="564nm_redu_O"+str(x),
                                                                                  ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax)
        except:
            print("Error "+str(x))
        else:
            print("Success "+str(x))
