# hantavirus-evolution
Sequence-based phylogenetic analysis and host prediction for hantavirus strains.

## Project Goal (we don't need all of these by midterms)
Build a pipeline that:

- **Downloads and cleans hantavirus sequences** - We grab all available 
hantavirus DNA sequences from a public government database, remove the ones 
with missing or bad data, and organize them so they're ready to analyze

- **Builds time calibrated phylogenetic trees** - We build a family tree of 
hantavirus strains showing how they're related to each other and how fast 
they've been changing over time

- **Identifies genomic sites under positive selection** - We find the specific 
spots in the virus's DNA that are evolving the fastest (not sure how possible this is, hantaviruses don't mutate that much compared to other virus families?)

- **Predicts reservoir host genus from sequence** - We train a machine learning 
model that can look at an unknown hantavirus sequence and predict which type 
of rodent it most likely came from (is this redundant?)

- **Flags fast evolving lineages for closer surveillance** - We automatically 
identify which hantavirus strains are changing unusually quickly compared to 
their relatives, so researchers know which ones to watch more closely

- **Visualizes all results in a live dashboard** - We build an interactive 
dashboard where you can see all our results, upload a new sequence, and get a 
prediction instantly without needing to run any code yourself (compiling speed would be really slow?)

## Team
- Aditri Mohanty
- Ana Pinkerton
- Srikritha Kosuri

## Project Structure
- Find research on the background of hantaviruses, how they spread, where they originate from, currrently documented strains, what researchers have done already, if a similar thing exists for other virus families, and project limiations (do this by monday?)
- afterwards we put it into the paper/presentation
