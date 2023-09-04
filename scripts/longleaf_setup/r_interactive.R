# Setting up vscode for R, longleaf, and interactive sessions

# referencing https://code.visualstudio.com/docs/languages/r
utils::install.packages('languageserver',
                 "/nas/longleaf/home/bgural/mambaforge/envs/cc_r_env/lib/R/library",
                 repos='http://cran.us.r-project.org' )
utils::install.packages('httpgd',
                 "/nas/longleaf/home/bgural/mambaforge/envs/cc_r_env/lib/R/library",
                 repos='http://cran.us.r-project.org' )

# installing to manage errors from the vscode R (Extension)
utils::install.packages('jsonlite'
                 "/nas/longleaf/home/bgural/mambaforge/envs/cc_r_env/lib/R/library",
                 repos='http://cran.us.r-project.org' )
