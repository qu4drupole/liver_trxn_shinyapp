oldRepos <- getOption("repos")
options(repos = c(oldRepos, 
                  BioCsoft = "https://bioconductor.org/packages/3.7/bioc",
                  BioCann = "https://bioconductor.org/packages/3.7/data/annotation",
                  BioCexp = "https://bioconductor.org/packages/3.7/data/experiment",
                  BioCworkflows = "https://bioconductor.org/packages/3.7/workflows"))
getOption("repos")



