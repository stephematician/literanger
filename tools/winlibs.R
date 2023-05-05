# Download eigen 3.4.0 headers

if (!file.exists("../winlibs/eigen3/Eigen/Eigen")) {

    download.file("https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip",
                  "eigen-3.4.0.zip",
                  quiet=T)
    dir.create("../winlibs", showWarnings=F)
    unzip("eigen-3.4.0.zip", exdir="../winlibs")
    file.rename('../winlibs/eigen-3.4.0', '../winlibs/eigen3')
    unlink("eigen-3.4.0.zip")

}

