# Download eigen 3.4.0 headers

if (!file.exists("../extlibs/eigen3/Eigen/Eigen")) {

    download.file("https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip",
                  "eigen-3.4.0.zip",
                  quiet=TRUE)
    dir.create("../extlibs", showWarnings=FALSE)
    unzip("eigen-3.4.0.zip", exdir="../extlibs")
    file.rename('../extlibs/eigen-3.4.0', '../extlibs/eigen3')
    unlink("eigen-3.4.0.zip")

}

