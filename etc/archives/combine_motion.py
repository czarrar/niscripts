#!/usr/bin/env python

import os, sys
import os.path as op

sys.path.append(join(environ.get("NISCRIPTS"), "include"))
from execute import Process    

def combine_motion(funcdir, runfile="good_run.txt", motionfile="motion.par", rundirs="run_*", 
                    outfile="motion_good.par"):
    fname = op.join(funcdir, rundirs, motionfile)
    motionfiles = glob(fname)
    if len(motionfiles) = 0:
        raise Exception("No motion files found for %s" % fname)
    motionfiles.sort()
    
    fname = op.join(funcdir, runfile)
    if not op.isfile(fname):
        raise Exception("No good runs file found for %s" % fname)
    f = file(fname, 'r')
    goodruns = f.readlines()
    f.close()
    if len(goodruns) == 0:
        print 'Empty good_run.txt in %s' % funcdir
        return 1
    goodruns = [ int(l.strip())-1 for l in goodruns if l.strip() ]
    
    motionfiles = np.array(motionfiles)[goodruns].tolist()
    
    outfile = op.join(funcdir, runfile)
    f = file(outfile, 'w')
    p = Process("cat %s" % " ".join(motionfiles), stdout=f, to_print=True)
    print tmp.stderr
    
    return 0

    
# 

argparse.ArgumentParser()