# slper -- SLiM Helper

For simulations that are compared to a single contemporary sample, it's likely
best to use TreeSeqeunces and compute summary statistics directly on those. For
simulations where you're tracking quantitative genetic simulation parameters
and allele frequencies every generation, `slper` is here to help.

This relies on using certain Eidos code blocks to record simulation in a
specific format, and then provides the Python code to process this output
efficiently.

## SLiM Code blocks

In the `initialize()` block, I have:

```
// a string for the file output name, including parameters
bname = ("simulation_name_" + N + "N_" + shift + "sd_" + rbp + 
         "rbp_" + alpha + "alpha_" +
         + mu + "mu_" + shifttime + "shifttime_" + seed);
defineConstant("basename", bname);
// quantgen variable output
defineConstant("stats_file", data_dir + '/' + basename + '_stats.tsv');

// metadata header for files
param_string = "#seed=" + seed + ";shift=" + shift + ";alpha=" + alpha + ";N=" + N + 
     	         ";rbp=" + rbp + ";region_length=" + region_length + 
               ";mu=" + mu + 
               ";Vs=" + Vs;

// write metadata header
writeFile(stats_file, param_string, append=F);
// write stats file column headers
writeFile(stats_file, "gen\tzbar\tzbar_nofix\tzvar\tgenic_var\tneut_het\tS", append=T);

// path to neutral frequencies
defineConstant("neutfreqs_file", data_dir + '/' + basename + '_neutfreqs.tsv');
// write neutral frequency metadata header
writeFile(neutfreqs_file, param_string, append=F);

```

Then, we data is recorded, we have code in the epoch blocks like the following
to record population-level statistics:

```
// mean population trait value
zbar = mean(inds.tagF); // mind your substitution policy!

// trait-affecting mutations
x = sim.mutationFrequencies(p1, sim.mutationsOfType(m1));
genic_var = 2*sum(alpha^2 * x * (1-x));

// neutral mutations
y = sim.mutationFrequencies(p1, sim.mutationsOfType(m2));
ssh = 2*sum(y * (1-y));

row = paste(c(sim.generation, zbar, pheno_mean, 
              var(inds.tagF), genic_var, ssh, length(y)), sep="\t");
writeFile(stats_file, row, append=T);

```

Then, to save frequencies, we need to following code blocks, which create a
ragged array output:


```
// the function to format frequencies
fmt_freqs='applyValue.id + ";" + applyValue.position + ";" + sim.mutationFrequencies(p1, applyValue);';

line = paste(sapply(sim.mutationsOfType(m2), fmt_freqs), sep="\t");
writeFile(neutfreqs_file, sim.generation + "\t" + line, append=T);
```

This all seems a bit verbose, mostly because currently we can't have Eidos
modules that are loaded in (otherwise I'd package more general versions of this
code there).


## Python Loading of SLiM Helper Results

The benefit is the way we can load this data easily into Python. Note that I
don't use the fully formatted filename as above, but just `results_*` for
brevity:

```
> stats = sf.parse_slim_stats('results_stats.tsv')
> stats.params
{'seed': '1',
	'shift': '1',
	'alpha': '0.01',
	'N': '1000',
	'rbp': '9e-09',
	'region_length': 50000000.0,
	'mu': '2e-07',
	'Vs': '1\n'}
> stats.stats.head()
gen zbar  zbar_nofix  zvar  genic_var neut_het  S
0 2 -0.00290  -0.00290  0.001994  0.002009  20.1966 16315
1 3 -0.00271  -0.00271  0.002931  0.002974  30.3371 20985
2 4 -0.00536  -0.00536  0.003923  0.003965  40.5791 24836
3 5 -0.00796  -0.00796  0.004793  0.004939  50.4836 27662
4 6 -0.00177  -0.00177  0.005581  0.005903  60.6699 30650
```

Then to load the frequency data:

```
> nfq = sf.parse_slim_ragged_freqs('results_neutfreqs.tsv')
> nfq.freqs.shape
(301, 42054823)
> nfq.params
{'seed': '1',
	'shift': '1',
	'alpha': '0.01',
	'N': '1000',
	'rbp': '9e-09',
	'region_length': 50000000.0,
	'mu': '2e-07',
	'Vs': '1\n'}
> nfq.samples  # timepoints saved
[1000,
 10001,
 ...]
```
