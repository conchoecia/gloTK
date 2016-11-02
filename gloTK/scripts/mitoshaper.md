This program:
  - assembles mitochondrial genomes using a seed-and-bait approach.
  - is useful in situations in which the reference and subject
    mitochondrial genomes may have different structures. This is a
    common phenomenon in phyla outside of vertebrates.
  - is not necessarily useful when mitochondrial structure is very
    similar, as one sees in the variations between human mitogenomes.

Assembly steps:
  1. sketch
    - Cleans up reads using Seqprep2
    - Downsamples the data
    - Performs a preliminary assembly with MIA
    - makes a model of mapped reads
  2. shape
    - Takes output model of mitoshaper-sketch and resamples data
    - determines seeds based on model
    - runs MitoBIM with each seed
  3. laminate
    - aligns all of the "plys" (contigs) from mitoshaper-shape
    - makes a kmer synteny matrix using the reference genome
    - finds the correct coordinates for the new assembly
  4. coat
    - Feeds the new assembly into MIA using unmerged, resampled reads
    - Generate skyline and synteny plots for the first iteration
  5. finish
    - Attempts to fix misassemblies
    - generates HTML reports for each assembly
