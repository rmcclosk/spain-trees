#!/usr/bin/env python

# approximate Bayesian computation (ABC) embedded in Markov chain Monte Carlo (MCMC)
# to estimate selection coefficients of alleles
# where alleles are defined as bins of g2p score
import csv
from scipy import optimize
from scipy import stats

# read in data
reader = csv.DictReader(open('00_raw.csv'))

data = {}
for row in reader:
    patid = row["Patient"]
    time = int(row["Time point"].split()[1])
    
    if not patid in data: data.update({patid:{}})
    if not time in data[patid]: data[patid].update({time:{}})
    
    try:
        fpr = float(row["G2P_FPR%"])
    except ValueError:
        continue

    if not fpr in data[patid][time]: data[patid][time].update({fpr:0.})
    data[patid][time][fpr] += int(row["Count"])

# normalize counts by totals
for patid in data.keys():
    for time in data[patid].keys():
        total_count = sum([c for c in data[patid][time].values()])
        for g2p in data[patid][time].keys():
            data[patid][time][g2p] = data[patid][time][g2p] / total_count

bins = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 25, 50, 100]
nbins = len(bins)-1    # each bin is defined by its left cutpoint

# collapse observed g2p FPR values into bins
bindata = {}
for patid in data.keys():
    bindata.update( { patid : {} } )
    for time in data[patid].keys():
        bindata[patid].update( { time : [0 for i in range(nbins)] } )
        for i in range(nbins):
            for fpr, freq in data[patid][time].items():
                if bins[i] <= fpr and fpr < bins[i+1]:
                    bindata[patid][time][i] += freq

# use cumulative distribution function of gamma as fitness function
#    shape > 1 and scale < 1 produces a bell-shaped curve that results in
#    a sigmoidal CDF.  use location parameter to shift around the midpoint
gamfun = stats.gamma

def fitness (x, params):
    shape, loc, scale, alpha = params
    res = alpha * gamfun.cdf(x, shape, loc, scale).tolist()
    return 1-res

def simulate (data, params):
    shape, loc, scale, alpha = params
    
    # fitness of different alleles
    fits = [ fitness(bins[i], params) for i in range(nbins) ]
    
    # return nested lists, one for each time point
    # use empty list to denote missing time point
    # cardinality of internal list is number of bins (n)
    sim = {}
    
    for patid in data.iterkeys():
        last_day = max(data[patid].keys())+1
        sim.update( {patid:[ [0 for i in range(nbins)] for d in range(last_day) ]} )
        
        # copy over initial allele frequencies - these don't change
        for i in range(nbins):
            sim[patid][0][i] = data[patid][0][i]
        
        
        for day in range(1, last_day):
            # mean fitness for previous allele frequency distribution
            meanfit = sum([ sim[patid][day-1][i] * fits[i] for i in range(nbins)])
            deltap = [ sim[patid][day-1][i] * (fits[i] - meanfit) / meanfit for i in range(nbins)]
            sim[patid][day] = [ sim[patid][day-1][i] + deltap[i] for i in range(nbins) ]
        
    return sim


def epsilon (data, sim):
    # compare simulation to observed data using least squares
    eps = 0.
    
    for patid in data.iterkeys():
        for time in data[patid].iterkeys():
            if time == 0:
                continue    # don't bother comparing initial allele frequencies
            # calculate squared distance
            eps += sum([ (data[patid][time][i]-sim[patid][time][i])**2 for i in range(nbins) ])
    
    return eps



# =================================================

def objfunc (x):
    shape, loc, scale, alpha = x
    if shape < 0:
        # bad parameter values, spike optimizer
        eps = 10000
    eps = epsilon(bindata, simulate(bindata, x))
    return eps

res = optimize.fmin(objfunc, (10, 0.05, 2, 0.4))
print('shape <- {}\nlocation <- {}\nscale <- {}\nalpha <- {}'.format(*res))
