# `sampleMass` output format

`sampleMass` outputs two files: a `<prefix>.massSamples` file and a `<prefix>.membership` file, both with the prefix specified in your base9.yaml.


## `<prefix>.massSamples`

Output format is two columns per row per star in the input photometry. The pairs of columns represent tuples of primary mass and secondary mass ratio. Each row represents a sample row in the `<prefix>.res` file.


## `<prefix.massSamples`

Output format is one column per row per star in the input photometry. Each column represents the cluster mebership likelihood of a given star. Each row represents a sample row in the `<prefix>.res` file.
