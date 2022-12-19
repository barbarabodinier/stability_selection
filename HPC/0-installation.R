setwd("../Package")

untar("fake_1.3.0.tar")
devtools::install("fake", upgrade = "always")

untar("sharp_1.2.1.tar.gz")
devtools::install("sharp", upgrade = "always")

devtools::install_version("glmnet", version = "4.1-1")

print(paste0("fake: ", packageVersion("fake")))
print(paste0("glmnet: ", packageVersion("glmnet")))
print(paste0("sharp: ", packageVersion("sharp")))
