# load everything required by dynamic topmodel
rm(list=ls())

# this loads common libraries
source("~/source/init.r")

src_dir <- "~/source/dynamictopmodel/dev"
lib_dir <- "~/source/lib"
source.dir(src_dir)
fns <- dir(src_dir, "\\.+r", full.names=T, recursive=T)
sapply(fns, source)
# shared source
source("~/source/lib/read_obs.r")
# p.e. calculations
source("~/source/lib/evap.r")

# support files
fns <- dir(file.path(src_dir, "dev", "*.r", recursive=T, full.names=T)







