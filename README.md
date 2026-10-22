# hantavirus-evolution
Sequence-based phylogenetic analysis and host prediction for hantavirus strains.

## Project Goal
Build a pipeline that:

- **Downloads and cleans hantavirus sequences** - We grab all available 
hantavirus DNA sequences from a public database, then clean and organize them so they're ready to analyze

- **Builds time calibrated phylogenetic trees** - We build a family tree of 
hantavirus strains showing how they're related to each other and how fast 
they've been changing over time

- **Identifies genomic sites under positive selection** - We find the specific 
spots in the virus's DNA that are evolving the fastest

- **Predicts reservoir host genus from sequence** - We train a machine learning 
model that can look at an unknown hantavirus sequence and predict which type 
of rodent it most likely came from

- **Flags fast evolving lineages for closer surveillance** - We automatically 
identify which hantavirus strains are changing unusually quickly compared to 
their relatives, so researchers know which ones to watch more closely

- **Visualizes all results in a live dashboard** - We build an interactive 
dashboard where you can see all our results, upload a new sequence, and get a 
prediction instantly without needing to run any code yourself

## Team
- Aditri Mohanty - Pipeline Coding, Slides
- Ana Pinkerton - Background Research, Phylogenetic tree, Slides
- Srikritha Kosuri - Background Research, History, Slides
