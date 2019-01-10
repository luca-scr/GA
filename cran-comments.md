## Test environments

* local macOS 10.12 Sierra, R 3.5.2
* rhub Debian Linux, R-release
* rhub Debian Linux, R-devel
* win-builder, R-release
* win-builder, R-devel

## R CMD check results

- local macOS 10.12 Sierra R 3.5.2
0 errors | 0 warnings | 0 note

- Debian Linux R-release
0 errors | 0 warnings | 1 note

* checking installed package size ... NOTE
  installed size is  7.3Mb
  sub-directories of 1Mb or more:
    doc    1.7Mb
    libs   5.1Mb
    
- Debian Linux R-devel
0 errors | 0 warnings | 0 note

- Windows R-release
0 errors | 0 warnings | 0 note

- Windows R-devel
0 errors | 0 warnings | 0 note


## Downstream dependencies

I have also run reverse dependency checking using the revdepcheck package. 
All packages that I was able to install did not show any issues related to GA.

