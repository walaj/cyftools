library(R.matlab)
library(rhdf5)

# download files
if (!file.exists(dest40 <- "~/projects/orion/matlab-CRC_Immune_MC-20210826.mat"))
    download.file("https://www.dropbox.com/s/va8xrq4buyarbz8/matlab-CRC_Immune_MC-20210826.mat?dl=1", dest40)
if (!file.exists(dest74 <- "~/projects/orion/matlab-CRCWSI_Tumor_MC-20210817.mat"))
    download.file("https://www.dropbox.com/s/qt60jn7k9rm4m9g/matlab-CRCWSI_Tumor_MC-20210817.mat?dl=1", dest74) 

                                        # read
contents <- h5ls(dest40)
data <- h5read(dest40, "dataTNPCRC_01")
data <- readMat(dest40)
