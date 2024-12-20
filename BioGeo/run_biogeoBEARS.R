# Load the necessary libraries
library(BioGeoBEARS)
library(ape)

# Step 1: Create a simple phylogenetic tree (Newick format)
# A simple tree with 4 species: A, B, C, D
test_tree <- read.tree(text = "((A,B), (C,D));")

# Step 2: Create simple geographic range data
# 4 species, 2 areas (1 represents presence, 0 represents absence)
test_geog <- data.frame(
  A = c(1, 0),  # Species A is in area 1
  B = c(1, 0),  # Species B is in area 1
  C = c(0, 1),  # Species C is in area 2
  D = c(0, 1)   # Species D is in area 2
)

# Write the geographic data to a file in PHYLIP format
# Ensure the first line is the number of species and areas
write.table(test_geog, "test_geog.data", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

# Open the file and prepend the number of species and areas
geog_file <- readLines("test_geog.data")
geog_file <- c("4 2", geog_file)  # Add the number of species and areas to the first line
writeLines(geog_file, "test_geog.data")  # Write it back to the file

# Check the contents of the file to ensure it's correct
cat("Contents of the file after writing:\n")
cat(readLines("test_geog.data"), sep = "\n")

# Step 3: Create the BioGeoBEARS run object
# Define the file paths for the tree and geographic data
tree_file <- "test_tree.newick"  # Save the tree in Newick format
write.tree(test_tree, file = tree_file)  # Write the tree to a file

# Create the BioGeoBEARS run object with necessary parameters
run_object <- list()
run_object$trfn <- tree_file  # Path to the tree file
run_object$geogfn <- "test_geog.data"  # Path to the geographic data file
run_object$max_range_size <- 2  # Maximum number of areas a species can occupy
run_object$use_detection_model <- FALSE  # Set detection model to FALSE if not needed

# Step 4: Run the DEC model (Dispersal-Extinction-Cladogenesis) using the bears_optim_run() function
test_res <- bears_optim_run(run_object)

# Step 5: View the results
print(test_res$outputs)

# Optional: Visualize the results on the tree
plot_BioGeoBEARS_results(
  results_object = test_res, 
  analysis_titletxt = "BioGeoBEARS DEC Model",
  addl_params = list("j"), 
  plotwhat = "text", 
  cornercoords_loc = NULL
)
