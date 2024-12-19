# Example Phylogeny

trfn = np(paste(addslash(extdata_dir), "Psychotria_5.2.newick", sep=""))
moref(trfn)
pdffn = "tree.pdf"
pdf(file=pdffn, width=9, height=9)
tr = read.tree(trfn)
tr

#plot(tr, cex=0.5)
#title("Example Psychotria phylogeny from Ree & Smith (2008)")
#axisPhylo() # plots timescale
#mtext("Millions of years ago (Ma)", side=1, line=2)

#dev.off()
#cmdstr = paste0("open ", pdffn)
#system(cmdstr)

# Example Geography

geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))
moref(geogfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = 4

#running DEC+J

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

BioGeoBEARS_run_object$num_cores_to_use = 1

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Psychotria_DEC+J_M0_unconstrained_v1.Rdata"

runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}

res = bears_optim_run(BioGeoBEARS_run_object)

#Plotting DEC+J

analysis_titletxt ="BioGeoBEARS DEC+J on Psychotria M0_unconstrained"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

