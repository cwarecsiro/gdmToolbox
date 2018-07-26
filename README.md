## Tools to fix the gdmEngine with...

### Install  (windows only)
`library(devtools)`  
`install_github('cwarecsiro/gdmToolbox', quick = TRUE, upgrade_dependencies = FALSE)`  
(Note... `quick = TRUE` will among other things bypass complaints about installing on x86 architecture,
while stopping `upgrade_dependencies` should prevent needless breaking of installed and loaded packages)