import numpy as np
import msprime
import twisst
import gzip
import sys

###
### MODEL -----------------------------------------------------------------------------------------
###

### SINGLE ORIGIN HYBRIDS MODEL ###
#Two parent populations split and then produce a single hybrid population, which later splits into two

def sim_hybrids_shared(pop_n, pop_Ne_P1, pop_Ne_P2, pop_Ne_H1, pop_Ne_H2, pop_Ne_OG, pop_Ne_Anc, pop_Ne_P12_Anc, pop_Ne_H12_Anc, pop_Ne_resize_P1, pop_Ne_resize_P2, pop_Ne_resize_P1_0, pop_Ne_resize_P2_0, prop_P1_H, t_hyb, t_hyb_split, t_parents, t_outgroup, t_resize, mig_P2P1_ancestral, mig_P2P1_recent, mig_P1P2_recent, l, r):
    demography = msprime.Demography()
    
    # be sure the first 5 pops are P1 P2 H1 H2 OG in this order!
    demography.add_population(name="P1", initial_size=pop_Ne_P1)
    demography.add_population(name="P2", initial_size=pop_Ne_P2)
    demography.add_population(name="H1", initial_size=pop_Ne_H1)
    demography.add_population(name="H2", initial_size=pop_Ne_H2)
    demography.add_population(name="OG", initial_size=pop_Ne_OG)
    demography.add_population(name="H12_Anc", initial_size=pop_Ne_H12_Anc)
    demography.add_population(name="Anc", initial_size=pop_Ne_Anc)
    demography.add_population(name="P12_Anc", initial_size=pop_Ne_P12_Anc)
    
    # Add recent, 2-way migration between ancestral parental populations
    demography.set_migration_rate(source="P2", dest="P1", rate=mig_P2P1_recent)
    demography.set_migration_rate(source="P1", dest="P2", rate=mig_P1P2_recent)
    
    # add split that produces two hybrid populations
    demography.add_population_split(time=t_hyb_split, derived=["H1", "H2"], ancestral="H12_Anc")
    
    # admixture event to produce hybrid ancestor population
    demography.add_admixture(time=t_hyb, derived="H12_Anc", ancestral=["P1","P2"], proportions=[prop_P1_H,1-prop_P1_H])
    
    # Reset ancestral migration from ancestral aq to ancestral pol populations
    demography.add_migration_rate_change(time=t_resize, source = "P2", dest = "P1", rate = mig_P2P1_ancestral)
    demography.add_migration_rate_change(time=t_resize, source = "P1", dest = "P2", rate = 0)
    
    # Resize parental populations
    demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P1, population="P1")
    demography.add_population_parameters_change(time=t_resize, initial_size=pop_Ne_resize_P2, population="P2")
    
    # Resize parental populations before admixture (backwards in time)
    demography.add_population_parameters_change(time=t_hyb, initial_size=pop_Ne_resize_P1_0, population="P1")
    demography.add_population_parameters_change(time=t_hyb, initial_size=pop_Ne_resize_P2_0, population="P2")
    
    # add population splits
    demography.add_population_split(time=t_parents, derived=["P1", "P2"], ancestral="P12_Anc")
    demography.add_population_split(time=t_outgroup, derived=["P12_Anc", "OG"], ancestral="Anc")
    
    # sort events
    demography.sort_events()
    
    # Simulate a tree sequence
    ts = msprime.sim_ancestry(samples={"P1":pop_n, "P2":pop_n, "H1":pop_n, "H2":pop_n, "OG":1},
                              demography=demography, ploidy = 2, sequence_length = l, recombination_rate = r)
    return(ts)


###
### PARAMETERS ------------------------------------------------------------------------------------
###

# LanW (H1) - LanR (H2) RIGHT

pop_Ne_Anc=500000/2
pop_Ne_P12_Anc=430000/2
pop_Ne_H12_Anc=10/2
pop_Ne_resize_P1=310000/2
pop_Ne_resize_P2=210000/2
pop_Ne_H1=89/2
pop_Ne_H2=79/2
pop_Ne_OG=200000/2
pop_Ne_resize_P1_0=52000/2
pop_Ne_resize_P2_0=280000/2
pop_Ne_P1=133/2
pop_Ne_P2=54/2
prop_P1_H=0.39
t_hyb=20
t_hyb_split=18
t_parents=225000
t_outgroup=2000000 # in theory, around 2e6 gens: 5Mya (Goropash. 2012) * 2.5 years / generation
t_resize=7500
mig_P2P1_ancestral=5.99e-6
mig_P2P1_recent=1.14e-5
mig_P1P2_recent=4.02e-6
pop_n=10
pop_n2=10
r=1e-6
n_blocks=100
block_length=1e4 # 1e4

###
### SIMULATE --------------------------------------------------------------------------------------
###

# run simulations to produce tree sequence objects

ts_shared_blocks = [sim_hybrids_shared(pop_n, pop_Ne_Anc=pop_Ne_Anc, pop_Ne_P12_Anc=pop_Ne_P12_Anc, pop_Ne_H12_Anc=pop_Ne_H12_Anc, pop_Ne_resize_P1=pop_Ne_resize_P1, pop_Ne_resize_P2=pop_Ne_resize_P2, pop_Ne_resize_P1_0=pop_Ne_resize_P1_0, pop_Ne_resize_P2_0=pop_Ne_resize_P2_0, pop_Ne_H1=pop_Ne_H1, pop_Ne_H2=pop_Ne_H2, pop_Ne_OG=pop_Ne_OG, pop_Ne_P1=pop_Ne_P1, pop_Ne_P2=pop_Ne_P2, prop_P1_H=prop_P1_H, t_hyb=t_hyb, t_hyb_split=t_hyb_split, t_parents=t_parents, t_outgroup=t_outgroup, t_resize=t_resize, mig_P2P1_ancestral=mig_P2P1_ancestral, mig_P2P1_recent=mig_P2P1_recent, mig_P1P2_recent=mig_P1P2_recent, l=block_length, r=r) for i in range(n_blocks)]


###
### VCF -------------------------------------------------------------------------------------------
###

# add mutations

ts_shared_blocks_mutated = [None]*n_blocks
for i in range(n_blocks):
    ts_shared_blocks_mutated[i] = msprime.sim_mutations(ts_shared_blocks[i], rate=3.5e-9)
    
# write VCF file

with gzip.open("../output.vcf.gz", "wt") as vcf_file:
    for i in range(n_blocks):
        ts_shared_blocks_mutated[i].write_vcf(vcf_file)


###
### TWISST ----------------------------------------------------------------------------------------
###

#for weighting we need to give twisst the names of each population, which are just numbers stored in the ts object by default
#note that all defined populations (including ancestral ones) are stored, but we only name those that we sampled from
parent1,parent2,hybrid1,hybrid2,outgroup = "0","1","2","3","4"

# h1
weights_shared_blocks1 = [None]*n_blocks
for i in range(n_blocks):
    print("block", i)
    weights_shared_blocks1[i] = twisst.weightTrees(ts_shared_blocks[i], treeFormat="ts",
                                               taxonNames=[parent1, parent2, hybrid1, outgroup],
                                               outgroup=outgroup, verbose=False)

twisst.summary(weights_shared_blocks1[0])

# h2
weights_shared_blocks2 = [None]*n_blocks
for i in range(n_blocks):
    print("block", i)
    weights_shared_blocks2[i] = twisst.weightTrees(ts_shared_blocks[i], treeFormat="ts",
                                               taxonNames=[parent1, parent2, hybrid2, outgroup],
                                               outgroup=outgroup, verbose=False)


#write to output files
with gzip.open("../sim1.tsv.gz", "wt") as weights_file:
    for i in range(n_blocks):
        twisst.writeWeights(weights_file, weights_shared_blocks1[i], include_topologies=True if i==0 else False, include_header=True if i==0 else False)

with gzip.open("../sim2.tsv.gz", "wt") as weights_file:
    for i in range(n_blocks):
        twisst.writeWeights(weights_file, weights_shared_blocks2[i], include_topologies=True if i==0 else False, include_header=True if i==0 else False)
