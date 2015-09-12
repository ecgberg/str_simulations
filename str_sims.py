'''
 ~~~~~ !!!!! ~~~~~ STR SIMULATION CODE v3 ~~~~~ !!!!! ~~~~~
 Takes in a population size, a number of generations to run for, mutational parameters, mutation effect size, and fitness parameters
 Initializes either a SNP or an STR population, based on the mutational parameters ([mu_SNP] for SNPs or [mu_STR, p, sigsq_g] for STRs)
 For each generation, pair chromosomes (no recombination for now), calculate fitness (attach to each chromosome), and output genotype/phenotype distributions
 Then WF sample *chromosomes* according to fitness (by keys) and then select chromosomes to mutate
 
 SNP chromosomes, with their associated genotypes, are stored as dictionaries; 
    the keys are the positions and the values are +-1 (the dosage)
 STR chromosomes are stored as values

 Information about a chromosome's allelic dosage, effect on gene expression, and fitness (paired with another chromosome)
 are calculated using numpy arrays

 It is from a numpy array, therefore, that we select a list of 2N chromosome names to survive (proportional to their fitness)
 These 2N chromosomes can then pick up mutations
'''
import numpy as np
import copy
import math
import argparse

# THIS IS FOR TESTING - COMMENT IT OUT WHEN ACTUALLY RUNNING
#global generation
#generation=100
#global writeout
#writeout=100

'''
 Takes a dictionary of chromosome dictionaries and the effect size of each mutation, calculates the gene expression
 of those chromosomes.
 Returns these values in a 2D numpy array (chromosome name = key of the chromosome of interest, expression)
 Column 1 of this array is the expression values ([:,1])
 Note that this function WOULD alter the underlying chromosomes dictionary - which is why I've done the deep copy
'''
def get_expression(chromosomes, effect_size):
    chrs=copy.deepcopy(chromosomes)
    for chrom, loci in chrs.iteritems():
        # If the genotypes are SNP genotypes, sum to get the total allelic dosage for that chromosome
        # If not, move on
        try:
            loci=sum(loci.values())
        except AttributeError:
            pass
        # Add the effect size to the total allelic dosage
        chrs[chrom]=loci*effect_size
    # Format into a 2D numpy array
    expressions=zip(chrs.keys(), chrs.values())
    exp=np.array(expressions)
    return(exp)


'''
 Calculate the fitness of one individual from an optimal value
 and a variance in fitness
 
 Note that this fitness is NOT normalized to be a probability distribution (max fitness will be 1)
'''
def get_fit(phenotype, opt, sigsq_f):
    lamb=float(phenotype-opt)
    exponent=(lamb ** 2)/(2*sigsq_f)
    
    prob_survival=(math.e ** (-1*exponent))
    return(prob_survival)


'''
 Vectorize the above function to calculate the fitness from one phenotype
 To apply this to a vector of phenotypes
'''
get_fitnesses=np.vectorize(get_fit, otypes=[np.float])



'''
 Takes a dictionary of chromosomes, an effect size for each mutation, and three files, open for writing
 The first will be for writing out the phenotype distribution of each chromosome
 That will have to suffice until I figure out how to write genotype distributions
 The second will write out the individual phenotype distribution
 The third will write out the distribution of individual fitnesses
'''
def get_phens_fits(chroms, b, o, s_f, f_p, f_f, f_c):
    exps=get_expression(chroms, b)
    
    # Sum pairs of chromosomes to get total gene expression for sets of individuals
    phens=np.array([sum(tup) for tup in zip(exps[:,1][0::2], exps[:,1][1::2])])
    
    # Pass list of phenotypes to calculate fitnesses
    # For chromosomes with less than a 1 in 100 chance of leaving a single offspring,
    # set the fitness to 0 (lethal)
    fits=get_fitnesses(phens, o, s_f)
    near_lethal = fits < 0.01*(len(chroms) ** -1)
    fits[near_lethal] = 0
    
    # Normalize fitnesses to sum to one
    fits=fits/(sum(fits))
    
    if np.isnan(sum(fits)):
        return
    
    # Write phenotype and fitness distributions to a file
    if generation % writeout == 0:
        np.savetxt(f_c, exps[:,1].reshape(1, exps[:,1].shape[0]), fmt='%.2f')
        np.savetxt(f_p, phens.reshape(1, phens.shape[0]), fmt='%.2f')
        np.savetxt(f_f, fits.reshape(1, phens.shape[0]), fmt='%.5f')
    
    # concatenate into a 2D numpy array (column 1 is chromosomal gene expression, column 2 is individual phenotype
    # column 3 is individual fitness (individual is the chromosome next to you
    # divide fitnesses by two to account for each chromosome having equal probability of survival
    mat=np.c_[exps, np.repeat(phens,2), np.repeat(fits,2)/2]
    return(mat)



'''
'''
def wf_sample(tmin1, pop_size, beta, opt, sigsq_f, f_genos, f_phens, f_fits, f_chroms, s):
    # Calculate fitnesses of the previous generation
    chrom_dat=get_phens_fits(chroms=tmin1, 
                             b=beta, 
                             o=opt, 
                             s_f=sigsq_f, 
                             f_p=f_phens, 
                             f_f=f_fits,
                             f_c=f_chroms)
    
    # Pick names of chromosomes to keep in proportion to their fitness
    try:
        if s:
            print 'SELECTION NOW'
            survivor_chrs=np.random.choice(chrom_dat[:,0], size=2*pop_size, p=chrom_dat[:,3])
        else:
            print 'NEUTRAL'
            survivor_chrs=np.random.choice(chrom_dat[:,0], size=2*pop_size)
    except TypeError:
        crash_string='POPULATION CRASH - Generation '+ str(generation)
        f_chroms.write(crash_string)
        f_genos.write(crash_string)
        f_phens.write(crash_string)
        f_fits.write(crash_string)
        return
    
    # Build a list of chromosomes appropriately
    chr_list=[tmin1[chrom].copy() for chrom in survivor_chrs]
    
    # Build a list of chromosome names
    new_pop=dict(zip(range(2*pop_size), chr_list))
    
    # Return the new list of chromosomes
    return(new_pop)



'''
 Takes in a dictionary of chromosomes and a list of mutation parameters (appropriate for SNP mutational model)
 mu is the probability that you add a new SNP to a chromosome - this may need some tweaking to account for
   the number of sites we think are on our chromosome of interest
 also note that mu, even though it is a single value, must be given as a list
 for purposes of distinguishing between SNP and STR chromosomes
 
 This function just updates the dictionary of chromosomes so it does not need to return anything
 However, it will write out the newly mutated genotypes to a file
'''
def mutate_SNP(chroms, mut_params, f_genos):
    mu, = mut_params
    
    # ~~~~~ #
    # First, draw all the mutations that will occur
    # ~~~~~ #
    
    num_mutations=np.random.poisson(mu*len(chroms))
    muts=np.random.choice([1,-1], size=num_mutations)
    
    # ~~~~~ #
    # Second, build new SNP mutations
    # ~~~~~ #
    
    snp_counts=dict()
    for chrom in chroms.values():
        for SNP in chrom.keys():
            try:
                snp_counts[SNP]+=1
            except KeyError:
                snp_counts[SNP]=1
    
    try:
        start_id=max(snp_counts.keys()) + 1
    except ValueError:
        start_id=0
    
    mut_ids=range(start_id, start_id+num_mutations)
    
    # ~~~~~ #
    # Third, assign those new SNP mutations to chromosomes
    # ~~~~~ #

    chr_to_mut=np.random.choice(chroms.keys(), size=num_mutations, replace=False)
    print('MUTATIONS: {0} occurred at positions {1} on chromosomes {2}').format(num_mutations, mut_ids, chr_to_mut)
    # Could maybe be done more efficiently if I need to
    for chrom, k, v in zip(chr_to_mut, mut_ids, muts):
        assert k not in chroms[chrom]
        chroms[chrom][k]=v
    
    # This may not be the most efficient way to handle writing genotypes to files if lots of SNPs pop up and go extinct...
    if generation % writeout == 0:
        genos_write=np.asarray(snp_counts.items()).T
        np.savetxt(f_genos, genos_write, fmt='%d')


'''
 Takes in a dictionary of chromosomes and a list of mutation parameters (appropriate for STR mutational models)
 mu is the probability that a mutation occurs in an STR
 p is the probability that that mutation changes the length by exactly 1 unit (equal probs + or -)
 sigsq_g is the variance in the geometric distribution from which a multi-step mutation would be drawn (still equal + or -)
 
 This function just updates the dictionary of chromosomes so it does not need to return anything
 
 Also will write out the length alleles of this generation * beta if we are supposed to
'''
def mutate_STR(chroms, mut_params, f_genos):
    mu, p, sigsq_g = mut_params
    
    # ~~~~~ #
    # First, draw all the mutations that will occur
    # ~~~~~ #
    
    num_mutations=np.random.poisson(mu*len(chroms))
    
    # Decide how many mutations will be single step and how many will be multi-step
    # Calculate the parameter for the geometric distribution
    # Draw the appropriate number of multi-steps from the appropriate geometric distribution
    num_multis=sum(np.random.uniform(size=num_mutations)<=p)
    g_param=(np.sqrt(4*sigsq_g+1)-1)/(2*sigsq_g)
    multis=np.random.geometric(p=g_param, size=num_multis)
    
    # Build the complete set of STR mutations
    # Shuffle so they're in random order (this may not be totally necessary
    muts=np.r_[multis, np.ones(num_mutations-num_multis)]
    np.random.shuffle(muts)
    
    # Select direction of mutation (one for each mutation that will occur)
    direction=np.random.choice([1,-1], size=num_mutations)
    
    # Add the direction of effect
    muts=muts*direction
    
    # ~~~~~ #
    # Second, draw chromosomes to add mutations to
    # ~~~~~ #
    chr_to_mut=np.random.choice(chroms.keys(), size=num_mutations, replace=True)

    # ~~~~~ #
    # Third, add the new mutations to the appropriate chromosomes
    # ~~~~~ #
    for chr_id, mut in zip(chr_to_mut, muts):
        chroms[chr_id]=chroms[chr_id]+mut
    
    if generation % writeout == 0:
        genos=np.array(chroms.values())
        np.savetxt(f_genos, genos.reshape((1, genos.shape[0])), fmt='%d')




'''
'''
def get_outfiles(opt, fvar, mut_params, effect_size, p):
    if len(mut_params)==1:
        t='SNP.'+str(effect_size)
        mut='mu_'+str(mut_params[0])
    elif len(mut_params)==3:
        t='STR.'+str(effect_size)
        mut='mu_'+str(mut_params[0])+'.'+str(mut_params[1])+'.'+str(mut_params[2])
    fits='fit_'+str(opt)+'.'+str(fvar)
    suffix=t+'_'+mut+'_'+fits
    fnames=[p+'genos_'+suffix, p+'phenos_'+suffix, p+'fits_'+suffix, p+'chroms_'+suffix]
    return(fnames)

'''
'''
def init_pop(size, t):
    chrom_ids=range(size)
    # This is DEFINITELY not the most memory efficient way to initialize a SNP population...
    if len(t)==1:
        chrom_genos=[{} for _ in xrange(size)]
    elif len(t)==3:
        chrom_genos=np.zeros(size, dtype=int)
    pop=dict(zip(chrom_ids, chrom_genos))
    return(pop)


'''
'''
def main(N, optimal_expression, fitness_variance, mutational_parameters, mutational_effect, generations, w, path, sel):
    fname_genotypes, fname_phenotypes, fname_fitnesses, fname_chromosomes, = get_outfiles(opt=optimal_expression, 
                                                                                          fvar=fitness_variance, 
                                                                                          mut_params=mutational_parameters, 
                                                                                          effect_size=mutational_effect,
                                                                                          p=path)
    f_c=open(fname_chromosomes, 'w')
    f_g=open(fname_genotypes, 'w')
    f_p=open(fname_phenotypes, 'w')
    f_f=open(fname_fitnesses, 'w')
    
    global generation
    generation=0
    global writeout
    writeout=w
    
    # !! ACTUALLY REALLY NEED TO DEAL WITH THIS NOWW !! #
    pop=init_pop(size=2*N, t=mutational_parameters)
    
    while generation < generations:
        print('Generation {0} of {1}:').format(generation, generations)
    
        # Add mutations at the appropriate rate
        # If you get a type error while trying to add mutations, this is bc pop is 'None'
        # Which would happen when the population crashes. So, break
        try:
            mutate_STR(pop, mutational_parameters, f_g)
        except TypeError:
            break
        except ValueError:
            pass
        
        try:
            mutate_SNP(pop, mutational_parameters, f_g)
        except TypeError:
            break
        except ValueError:
            pass
        
        # Wright-Fisher sample
        pop=wf_sample(tmin1=pop, 
                        pop_size=N, 
                        beta=mutational_effect, 
                        opt=optimal_expression, 
                        sigsq_f=fitness_variance, 
                        f_genos=f_g, 
                        f_phens=f_p, 
                        f_fits=f_f,
                        f_chroms=f_c,
                        s=sel)
        
        generation=generation+1
    
    f_g.close()
    f_p.close()
    f_f.close()
    

parser=argparse.ArgumentParser(description='Acquire input variables')
parser.add_argument('-m',
                    type=float,
                    nargs='+',
                    help='Mutational parameters - for SNPS, [chrom mutation rate], for STRs, [chrom mutation rate, prob single step mutation, var in multi-mut step size')
parser.add_argument('-b',
                    type=float,
                    help='Linear effect size of each mutation')
parser.add_argument('-N',
                    type=int,
                    default=1000,
                    help='Number of individuals to include in each generation - default is 1000')
parser.add_argument('-o',
                    type=int,
                    default=0,
                    help='Optimal expression level - default is 0')
parser.add_argument('-f',
                    type=float,
                    default=1,
                    help='Variance in fitness around optimal - default is 1')
parser.add_argument('-g',
                    type=int,
                    default=0,
                    help='Number of generations to run simulation - default is 10N')
parser.add_argument('-w',
                    type=int,
                    default=100,
                    help='Interval between which to write out population information - default is 100')
parser.add_argument('-p', 
                    type=str, 
                    default='/Users/eglassbe/Dropbox/Pritchard_Lab/str_simulations/str_res/',
                    help='Path to write results files to - default is Pritchard Lab Dropbox str_res/')
parser.add_argument('--neutral',
                    dest='s',
                    action='store_false',
                    help='Turn off selection')
parser.add_argument('--no-neutral',
                    dest='s',
                    action='store_true',
                    help='Keep selection on')
parser.set_defaults(s=True)
                    
args=parser.parse_args()
print(args)

if args.g==0:
    args.g=10*args.N

if __name__ == '__main__':
    main(N=args.N, 
         optimal_expression=args.o, 
         fitness_variance=args.f, 
         mutational_parameters=args.m, 
         mutational_effect=args.b, 
         generations=args.g, 
         w=args.w,
         path=args.p,
         sel=args.s)
